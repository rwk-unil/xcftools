/*******************************************************************************
 * Copyright (C) 2024 Rick Wertenbroek
 * Copyright (C) 2022-2023 Olivier Delaneau
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

#include <modes/binary2sapphire.h>
#include <utils/xcf.h>

#include <containers/bitvector.h>
#include <objects/sparse_genotype.h>

#include <cmath>

using namespace std;

binary2sapphire::binary2sapphire(string _region, int _nthreads, float maf_threshold,
				 size_t fifo_size, bool pp_from_maf, bool pp_from_af) :
	region(_region), nthreads(_nthreads),
	FIFO_SIZE(fifo_size),
	PP_THRESHOLD(0.99),
	MAF_THRESHOLD(maf_threshold),
	pp_arr(NULL),
	pp_arr_size(0),
	start_id(0), stop_id(-1),
	line_counter(0),
	print_counter(0),
	pred(PP_THRESHOLD),
	progress(0),
	pp_from_maf(pp_from_maf),
	pp_from_af(pp_from_af) {
	if (pp_from_maf) {
		std::cout << "The PP score will be generated from MAF" << std::endl;
	}
	if (pp_from_af) {
		std::cout << "The PP score will be generated from AF" << std::endl;
	}
}

binary2sapphire::~binary2sapphire() {
}

void binary2sapphire::convert(string finput, string foutput) {
	tac.clock();

	vrb.title("Extracting from XCF to SAPPHIRE Binary");
	if (region.empty()) vrb.bullet("Region        : All");
	else vrb.bullet("Region        : " + stb.str(region));

	//Opening XCF reader for input
	xcf_reader XR(region, nthreads);
	int32_t idx_file = XR.addFile(finput);

	//Get file type
	int32_t type = XR.typeFile(idx_file);
	if (type != FILE_BINARY) vrb.error("[" + finput + "] is not a XCF file");

	//Get sample IDs
	vector < string > samples;
	int32_t nsamples = XR.getSamples(idx_file, samples);
	vrb.bullet("#samples = " + stb.str(nsamples));

	//Buffer for input/output
	int32_t * input_buffer = (int32_t*)malloc(2 * nsamples * sizeof(int32_t));
	int32_t * output_buffer = (int32_t*)malloc(2 * nsamples * sizeof(int32_t));

	//Buffer for binary data
	bitvector binary_buffer = bitvector(2 * nsamples);

	//Prepare parameters for SAPPHIRE binary
	number_of_het_sites.clear();
	number_of_low_pp_sites.clear();
	number_of_non_snp.clear();
	number_of_het_sites.resize(nsamples, 0);
	number_of_low_pp_sites.resize(nsamples, 0);
	number_of_snp_low_pp_sites.resize(nsamples, 0);
	number_of_non_snp.resize(nsamples, 0);

	/* Handle garbage input and limit to number of samples */
	if (start_id > nsamples) {
		start_id = nsamples;
	}
	if (stop_id > nsamples) {
		stop_id = nsamples;
	}
	if (stop_id < start_id) {
		stop_id = start_id;
	}
	std::cout << "Start ID : " << start_id << " Stop ID : " << stop_id << std::endl;
	fifos.resize(stop_id-start_id, GenericKeepFifo<HetInfo, PPPred>(FIFO_SIZE, PPPred(PP_THRESHOLD)));

	//Proceed with conversion
	uint32_t n_lines = 0;
	while (XR.nextRecord()) {

		//Get type of record
		type = XR.typeRecord(idx_file);

		//Convert from BCF; copy the data over
		if (type == RECORD_BCFVCF_GENOTYPE) {
			XR.readRecord(idx_file, reinterpret_cast< char** > (&input_buffer));
			memcpy(output_buffer, input_buffer, 2 * nsamples * sizeof(int32_t));
		}

		//Convert from binary genotypes
		else if (type == RECORD_BINARY_GENOTYPE) {
			int n = XR.readRecord(idx_file, reinterpret_cast< char** > (&binary_buffer.bytes));
			for(uint32_t i = 0 ; i < nsamples ; i++) {
				bool a0 = binary_buffer.get(2*i+0);
				bool a1 = binary_buffer.get(2*i+1);
				if (a0 == true && a1 == false) {
					output_buffer[2*i+0] = bcf_gt_missing;
					output_buffer[2*i+1] = bcf_gt_missing;
				} else {
					output_buffer[2*i+0] = bcf_gt_unphased(a0);
					output_buffer[2*i+1] = bcf_gt_unphased(a1);
				}
			}
		}

		//Convert from binary haplotypes
		else if (type == RECORD_BINARY_HAPLOTYPE) {
			XR.readRecord(idx_file, reinterpret_cast< char** > (&binary_buffer.bytes));
			for(uint32_t i = 0 ; i < nsamples ; i++) {
				bool a0 = binary_buffer.get(2*i+0);
				bool a1 = binary_buffer.get(2*i+1);
				output_buffer[2*i+0] = bcf_gt_phased(a0);
				output_buffer[2*i+1] = bcf_gt_phased(a1);
			}
		}

		//Convert from sparse genotypes
		else if (type == RECORD_SPARSE_GENOTYPE) {
			int32_t n_elements = XR.readRecord(idx_file, reinterpret_cast< char** > (&input_buffer)) / sizeof(int32_t);
			//Set all genotypes as major
			bool major = (XR.getAF()>0.5f);
			std::fill(output_buffer, output_buffer+2*nsamples, bcf_gt_unphased(major));
			//Loop over sparse genotypes
			for(uint32_t r = 0 ; r < n_elements ; r++) {
				sparse_genotype rg;
				rg.set(input_buffer[r]);
				if (rg.mis) {
					output_buffer[2*rg.idx+0] = bcf_gt_missing;
					output_buffer[2*rg.idx+1] = bcf_gt_missing;
				} else {
					output_buffer[2*rg.idx+0] = bcf_gt_unphased(rg.al0);
					output_buffer[2*rg.idx+1] = bcf_gt_unphased(rg.al1);
				}
			}
		}

		//Convert from sparse haplotypes
		else if (type == RECORD_SPARSE_HAPLOTYPE) {
			int32_t n_elements = XR.readRecord(idx_file, reinterpret_cast< char** > (&input_buffer)) / sizeof(int32_t);
			//Set all genotypes as major
			bool major = (XR.getAF()>0.5f);
			std::fill(output_buffer, output_buffer+2*nsamples, bcf_gt_phased(major));
			//Loop over sparse genotypes
			for(uint32_t r = 0 ; r < n_elements ; r++) output_buffer[input_buffer[r]] = bcf_gt_phased(!major);
		}

		//Unknown record type
		else vrb.bullet("Unrecognized record type [" + stb.str(type) + "] at " + XR.chr + ":" + stb.str(XR.pos));

		//SAPPHIRE Extraction
		bool has_pp = false;
		bool non_snp = false;
		if (XR.ref.size() > 1 || XR.alt.size() > 1) {
			non_snp = true;
		}

		// We need AC (Allele Count) and AN (Allele Number) for PP from MAF
		float synthetic_pp = 0.0;

		if (pp_from_maf) {
			float maf = float(XR.getAC())/XR.getAN();
			synthetic_pp = (maf > MAF_THRESHOLD) ? NAN : 0.5 + maf / 2.0;
			//std::cout << "AC : " << XR.getAC() << "\tAN : " << XR.getAN() << "\tsynth PP " << synthetic_pp << std::endl;
		}

		if (pp_from_af) {
			float AF = XR.getAF();
			if (AF > 0.5) {
				// Complement to get the minor allele frequency
				AF = 1.0 - AF;
			}
			synthetic_pp = (AF > MAF_THRESHOLD) ? NAN : 0.5 + AF / 2.0;
			//std::cout << "AF: " << AF << " Synth PP : " << synthetic_pp << std::endl;
		}

		// Extract heterozygous sites and PP
		for (size_t i = start_id; i < stop_id; ++i) {
			int encoded_a0 = output_buffer[i*PLOIDY_2];
			int encoded_a1 = output_buffer[i*PLOIDY_2+1];
			int a0 = bcf_gt_allele(encoded_a0);
			int a1 = bcf_gt_allele(encoded_a1);

			if (a0 != a1) {
				number_of_het_sites[i]++;
				float pp = NAN; // Assume perfect phasing if there is no PP field
				if (has_pp) {
					pp = pp_arr[i];
				} else if (pp_from_maf || pp_from_af) {
					pp = synthetic_pp;
				}

				if (XR.getAC() == 1 && pp >= PP_THRESHOLD) {
					/* Edge case for old version of SHAPEIT5 that would score
					singletons phased with only one of their parents with
					a PP of 1.0, In about 95% of cases the singleton comes
					from the other parent if not observed in the known parent
					but it could also be a de novo mutation on the known
					parent haplotype, these cases should not be scored 1.0 */
					pp = 0.97; /* arbitrary value */
				}

				HetInfo hi(line_counter, encoded_a0, encoded_a1, pp);

				if (pred(hi)) {
					number_of_low_pp_sites[i]++;
					if (!non_snp) {
						number_of_snp_low_pp_sites[i]++;
					}
				}
				if (non_snp) {
					number_of_non_snp[i]++;
				}

				fifos[i-start_id].insert(hi);
			}
		}

		line_counter++;
		if (progress) {
			if (++print_counter == progress) {
				print_counter = 0;
				printf("\033[A\033[2K");
				std::cout << "Handled " << line_counter << " VCF entries (lines)" << std::endl;
			}
		}

		//Counting
		n_lines++;

		//Verbose
		if (n_lines % 10000 == 0) vrb.bullet("Number of XCF records processed: N = " + stb.str(n_lines));
	}

	vrb.bullet("Number of XCF records processed: N = " + stb.str(n_lines));

	//Free
	free(input_buffer);
	free(output_buffer);

	//Close files
	XR.close();

	// Write SAPPHIRE Binary file
	finalize();
	show_info();
	write_to_file(foutput);
}

