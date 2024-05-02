/*******************************************************************************
 * Copyright (C) 2024 Rick Wertenbroek
 *
 * Permission is hereby granted, free of charge, to any person obtaining a copy
 * of this software and associated documentation files (the "Software"), to deal
 * in the Software without restriction, including without limitation the rights
 * to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
 * copies of the Software, and to permit persons to whom the Software is
 * furnished to do so, subject to the following conditions:
 *
 * The above copyright notice and this permission notice shall be included in
 * all copies or substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
 * AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
 * OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
 * SOFTWARE.
 ******************************************************************************/

#include "../versions/versions.h"
#include "modes/binarysapphire2binary.h"
#include "utils/xcf.h"

#include "containers/bitvector.h"
#include "objects/sparse_genotype.h"

#include "utils/sapphire/fs.hpp"
#include "utils/sapphire/het_info.hpp"
#include "utils/sapphire/het_info_loader.hpp"

class VCFLineWork {
public:
	VCFLineWork() : vcf_line_num(0) {}
	VCFLineWork(size_t vcf_line_num) : vcf_line_num(vcf_line_num) {}
	VCFLineWork(HetInfo& hi, size_t id) : vcf_line_num(hi.vcf_line) {
		updated_data[id] = hi;
	}

	size_t vcf_line_num;
	std::map<size_t, HetInfo> updated_data;
};

void insert_in_work(std::map<size_t, VCFLineWork>& work, HetInfo& hi, size_t id) {
	if (auto work_line = work.find(hi.vcf_line); work_line != work.end()) {
		work_line->second.updated_data[id] = hi;
	} else {
		work[hi.vcf_line] = VCFLineWork(hi, id);
	}
}

void fill_work_from_himm(std::map<size_t, VCFLineWork>& work, HetInfoMemoryMap& himm) {
	// For all samples
	for (size_t i = 0; i < himm.num_samples; ++i) {
		// Get Het variants
		std::vector<HetInfo> his;
		himm.fill_het_info(his, i);

		// For all het variants
		for (auto& hi : his) {
		// If it has been rephased (>1.0)
			if (!std::isnan(hi.pp) && hi.pp > 1.0) {
				// Insert in work
				insert_in_work(work, hi, i);
			}
		}
	}
}

binarysapphire2binary::binarysapphire2binary(std::string _region, float _minmaf, int _nthreads, int _mode, bool _drop_info)
{
	mode = _mode;
	nthreads = _nthreads;
	region = _region;
	minmaf = _minmaf;
	drop_info = _drop_info;
}

binarysapphire2binary::~binarysapphire2binary()
{
}

int32_t binarysapphire2binary::parse_genotypes(xcf_reader& XR, const uint32_t idx_file)
{
	//Get type of record
	const int32_t type = XR.typeRecord(idx_file);
	int32_t n_elements = XR.ind_names[idx_file].size();
	if (type == RECORD_BCFVCF_GENOTYPE) {
		vrb.error("BCF/VCF record in binary2binary mode !");
		//XR.readRecord(idx_file, reinterpret_cast< char** > (&sparse_int_buf));
	}
	else if (type == RECORD_BINARY_GENOTYPE) {
		XR.readRecord(idx_file, reinterpret_cast< char** > (&binary_bit_buf.bytes));
	}
	else if (type == RECORD_BINARY_HAPLOTYPE) {
		XR.readRecord(idx_file, reinterpret_cast< char** > (&binary_bit_buf.bytes));
	}
	else if (type == RECORD_SPARSE_GENOTYPE) {
		n_elements = XR.readRecord(idx_file, reinterpret_cast< char** > (&sparse_int_buf)) / sizeof(int32_t);
	}
	else if (type == RECORD_SPARSE_HAPLOTYPE) {
		n_elements = XR.readRecord(idx_file, reinterpret_cast< char** > (&sparse_int_buf)) / sizeof(int32_t);
	}
	else vrb.bullet("Unrecognized record type [" + stb.str(type) + "] at " + XR.chr + ":" + stb.str(XR.pos));

	return n_elements;
}

void binarysapphire2binary::convert(std::string finput, std::string foutput, std::string fsapphire)
{
	tac.clock();
	switch (mode)
	{
		case CONV_BCF_BG: vrb.title("Converting from XCF+SAPPHIRE to XCF [Binary/Genotype]"); break;
		case CONV_BCF_BH: vrb.title("Converting from XCF+SAPPHIRE to XCF [Binary/Haplotype]"); break;
		case CONV_BCF_SG: vrb.title("Converting from XCF+SAPPHIRE to XCF [Sparse/Genotype]"); break;
		case CONV_BCF_SH: vrb.title("Converting from XCF+SAPPHIRE to XCF [Sparse/Haplotype]"); break;
	}
	vrb.bullet("SAPPHIRE file : " + stb.str(fsapphire));

	if (region.empty()) {
		vrb.bullet("Region        : All");
	} else {
		vrb.bullet("Region        : " + stb.str(region));
		vrb.error("SAPPHIRE update does not suppport the region option");
	}

	if (mode == CONV_BCF_SG || mode == CONV_BCF_SH) vrb.bullet("Min MAF       : " + stb.str(minmaf));

	xcf_reader XR(1);
	const uint32_t idx_file = XR.addFile(finput);
	const int32_t typef = XR.typeFile(idx_file);
	if (typef != FILE_BINARY) vrb.error("[" + finput + "] is not a XCF file");
	uint32_t nsamples_input = XR.ind_names[idx_file].size();
	xcf_writer XW(foutput, false, nthreads);
	bcf1_t* rec = XW.hts_record;
	size_t line_counter = 0;
	size_t errors = 0;
	size_t updated_gts = 0;

	vrb.bullet("Generating workload from SAPPHIRE file...");

	HetInfoMemoryMap himm(fsapphire);
	std::map<size_t, VCFLineWork> work;
	fill_work_from_himm(work, himm);

	if (drop_info) XW.writeHeader(XR.sync_reader->readers[0].header, XR.ind_names[idx_file], std::string("XCFtools ") + std::string(XCFTLS_VERSION));
	else XW.writeHeaderClone(XR.sync_reader->readers[0].header,XR.ind_names[idx_file], std::string("XCFtools ") + std::string(XCFTLS_VERSION));

	binary_bit_buf.allocate(2 * nsamples_input);
	sparse_int_buf.resize(2 * nsamples_input,0);

	uint32_t n_lines_rare = 0, n_lines_comm = 0;

	while (XR.nextRecord())
	{
		//Is that a rare variant?
		float af =  XR.getAF();
		float maf = std::min(af, 1.0f-af);
		bool minor = (af < 0.5f);
		bool rare = (maf < minmaf);

		if (drop_info)
			XW.writeInfo(XR.chr, XR.pos, XR.ref, XR.alt, XR.rsid, XR.getAC(), XR.getAN());
		else
			XW.hts_record = XR.sync_lines[0];

		int32_t n_elements = parse_genotypes(XR,idx_file);
		const int32_t type = XR.typeRecord(idx_file);

		// Apply sapphire changes
		if (auto work_line = work.find(line_counter); work_line != work.end()) {
			if (work_line->second.vcf_line_num != line_counter) {
				vrb.error("Work line is different from line counter !");
				vrb.error("Something went wrong in the machinery");
				errors++;
			} else {
				for (auto todo : work_line->second.updated_data) {
					auto idx = todo.first;
					auto minor_idx0 = 2*idx;
					auto minor_idx1 = 2*idx+1;
	
					// Update GT
					if (type==RECORD_BINARY_GENOTYPE) {
						vrb.error("Binary genotype is for unphased data");
						errors++;
						// No need to update anything since het is always "01" ("10" is missing)
					} else if (type==RECORD_BINARY_HAPLOTYPE) {
						// Only update if rephased, no need otherwise
						if (binary_bit_buf.get(minor_idx0) != todo.second.a0) {
							updated_gts++;
							binary_bit_buf.set(minor_idx0, todo.second.a0);
							binary_bit_buf.set(minor_idx1, todo.second.a1);
						}
					} else if (type==RECORD_SPARSE_GENOTYPE) {
						vrb.error("Binary genotype is for unphased data");
						errors++;
						// No need to update anything since het is always "01" ("10" is missing)
					} else if (type==RECORD_SPARSE_HAPLOTYPE) {
						// We don't know if the original file has minor idx0 or minor idx1 set, so search which one
						// This may seem slow, but it is sparse because there are not many elements, so this find is quite fast
						auto it0 = std::find(sparse_int_buf.begin(), sparse_int_buf.end(), minor_idx0);
						auto it1 = std::find(sparse_int_buf.begin(), sparse_int_buf.end(), minor_idx1);

						if ((it0 != std::end(sparse_int_buf)) && (it1 != std::end(sparse_int_buf))) {
							vrb.error("Sample is hom alt, cannot rephase)");
							errors++;
						}
						if ((it0 == std::end(sparse_int_buf)) && (it1 == std::end(sparse_int_buf))) {
							vrb.error("Sample is hom ref, cannot rephase");
							errors++;
						} else {
							// Choose the correct iterator
							auto it = (it0 == std::end(sparse_int_buf)) ? it1 : it0;
							// Set the index of a0 if a0, else a1
							*it = (todo.second.a0 ? minor_idx0 : minor_idx1);
						}
					}
				}
			}
		}

		//Write record
		if (mode == CONV_BCF_SG && rare)
		{
			if (type==RECORD_SPARSE_GENOTYPE)
				XW.writeRecord(RECORD_SPARSE_GENOTYPE, reinterpret_cast<char*>(sparse_int_buf.data()), n_elements * sizeof(int32_t));
			else if (type==RECORD_BINARY_GENOTYPE)
			{
				//conversion: BINARY gen -> sparse
				n_elements=0;
				for(uint32_t i = 0 ; i < nsamples_input ; i++)
				{
					const bool a0 = binary_bit_buf.get(2*i+0);
					const bool a1 = binary_bit_buf.get(2*i+1);
					sparse_int_buf[n_elements++] = sparse_genotype(i, (a0!=a1), (a0 && !a1), a0, a1, 0).get();
				}
				XW.writeRecord(RECORD_SPARSE_GENOTYPE, reinterpret_cast<char*>(sparse_int_buf.data()), n_elements * sizeof(int32_t));
			}
			else vrb.error("Converting non-genotype type to genotype type!");
		}
		else if (mode == CONV_BCF_SH && rare)
		{
			if (type==RECORD_SPARSE_HAPLOTYPE)
				XW.writeRecord(RECORD_SPARSE_HAPLOTYPE, reinterpret_cast<char*>(sparse_int_buf.data()), n_elements * sizeof(int32_t));
			else if (type==RECORD_BINARY_HAPLOTYPE)
			{
				//conversion: BINARY hap -> sparse
				n_elements=0;
				for (size_t i = 0; i < 2 * nsamples_input; ++i)
				{
					if (binary_bit_buf.get(i) == true)
						sparse_int_buf[n_elements++]=i;
				}
				XW.writeRecord(RECORD_SPARSE_HAPLOTYPE, reinterpret_cast<char*>(sparse_int_buf.data()), n_elements * sizeof(int32_t));
			}
			else vrb.error("Converting non-haplotype type to haplotype type!");
		}
		else if (mode == CONV_BCF_SG || mode == CONV_BCF_BG)
		{
			if (type==RECORD_BINARY_GENOTYPE)
				XW.writeRecord(RECORD_BINARY_GENOTYPE, binary_bit_buf.bytes, binary_bit_buf.n_bytes);
			else if (type==RECORD_SPARSE_GENOTYPE)
			{
				binary_bit_buf.set(false);
				for (auto gt : sparse_int_buf)
				{
					sparse_genotype rg;
					rg.set(gt);
					if (rg.mis) {
						binary_bit_buf.set(2*rg.idx+0, true);
						//binary_bit_buf.set(2*rg.idx+1) = false;
					} else {
						binary_bit_buf.set(2*rg.idx+0,rg.al0);
						binary_bit_buf.set(2*rg.idx+1,rg.al1);
					}
				}
				XW.writeRecord(RECORD_BINARY_GENOTYPE, binary_bit_buf.bytes, binary_bit_buf.n_bytes);
			}
			else vrb.error("Converting non-genotype type to genotype type!");
		}
		else
		{
			if (type==RECORD_BINARY_HAPLOTYPE)
				XW.writeRecord(RECORD_BINARY_HAPLOTYPE, binary_bit_buf.bytes, binary_bit_buf.n_bytes);
			else if (type==RECORD_SPARSE_HAPLOTYPE)
			{
				//conversion: SPARSE hap -> binary
				binary_bit_buf.set(!minor);
				for (auto index : sparse_int_buf)
					binary_bit_buf.set(index,minor);
				XW.writeRecord(RECORD_BINARY_GENOTYPE, binary_bit_buf.bytes, binary_bit_buf.n_bytes);
			}
			else vrb.error("Converting non-haplotype type to haplotype type!");
		}
		//Line counting
		n_lines_comm += !rare || mode == CONV_BCF_BG || mode == CONV_BCF_BH;
		n_lines_rare += rare && (mode == CONV_BCF_SG || mode == CONV_BCF_SH);
		line_counter++;

		//Verbose
		if ((n_lines_comm+n_lines_rare) % 10000 == 0) {
			if (mode == CONV_BCF_BG || mode == CONV_BCF_BH) vrb.bullet("Number of BCF records processed: N=" + stb.str(n_lines_comm));
			else vrb.bullet("Number of BCF records processed: Nc=" + stb.str(n_lines_comm) + "/ Nr=" + stb.str(n_lines_rare));
		}
	}

	if (mode == CONV_BCF_BG || mode == CONV_BCF_BH) vrb.bullet("Number of records processed: N=" + stb.str(n_lines_comm));
	else vrb.bullet("Number of records processed: Nc=" + stb.str(n_lines_comm) + "/ Nr=" + stb.str(n_lines_rare));

	if (!drop_info) XW.hts_record = rec;

	XW.close();//always close XW first? important for multithreading if set
	XR.close();
}
