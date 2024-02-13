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

#ifndef _BINARY2SAPPHIRE_H
#define _BINARY2SAPPHIRE_H

#include <utils/otools.h>

#ifndef __HET_INFO_HPP__
#define __HET_INFO_HPP__

#include <string>
#include "vcf.h"

#ifndef __FS_HPP__
#define __FS_HPP__

#include <libgen.h> // Has dirname() / basename()
#include <string>
#include <cstring> // for strdup
#include <fstream>
#include <unistd.h>
/*
The functions dirname() and basename() break a null-terminated pathname string
into directory and filename components. In the usual case, dirname() returns the
string up to, but not including, the final '/', and basename() returns the
component following the final '/'. Trailing '/' characters are not counted as
part of the pathname.
*/

#if __cplusplus < 201703L
#define MAYBE_UNUSED
#warning "Warnings about unused functions are expected"
#else
#define MAYBE_UNUSED [[maybe_unused]]
#endif

/// @todo move this to an object file
namespace
{
	class NamedFileStream {
	public:
		NamedFileStream(std::string filename) : stream(filename, stream.binary | stream.out | stream.trunc), filename(filename) {}
		std::fstream stream;
		std::string filename;
	};

	MAYBE_UNUSED NamedFileStream get_temporary_file(int* file_desc, const char* name_template /** @todo flags */) {
		char *tmpname = strdup(name_template);
		int fd = mkstemp(tmpname);
		std::string filename(tmpname);
		free(tmpname);
		if (fd != -1) {
			NamedFileStream nfs(filename); // Opens as stream
			if (file_desc) {
				*file_desc = fd;
			} else {
				// Close the file descriptor if not used (file is open as stream)
				close(fd);
			}
			if (!nfs.stream.is_open()) {
				std::cerr << "Could not open temporary file stream to : " << filename << std::endl;
				throw "Could not get temporary file";
			}
			return nfs;
		} else {
			throw "Could not get temporary file";
		}

		return NamedFileStream(filename);
	}

	MAYBE_UNUSED NamedFileStream get_temporary_file(int* file_desc = nullptr /** @todo flags */) {
		return get_temporary_file(file_desc, "/tmp/tmpfileXXXXXX");
	}
} // namespace

#if __cplusplus >= 201703L
	#include <filesystem>
	namespace fs = std::filesystem;
	using fs::remove;
	// This can be used instead of basename() (more portable, if C++17 is available)
	// fs::path( "/foo/bar.txt" ).filename() => "bar.txt"
#else
	#include <stdio.h> // Has remove()
	#include <sys/stat.h>
	namespace fs {
		inline size_t file_size(const std::string& filename) {
			struct stat st;
			if (stat(filename.c_str(), &st) < 0) {
				std::cerr << "Size of file : " << filename << " could not be determined !" << std::endl;
				return 0;
			} else {
				return st.st_size;
			}
		}

		inline bool exists(const std::string& filename) {
			struct stat st;
			return stat(filename.c_str(), &st) == 0;
		}
	};
#endif

#undef MAYBE_UNUSED

#endif /* __FS_HPP__ */

class HetInfo {
public:
	HetInfo(std::ifstream& ifs) {
		from_stream(ifs);
	}
	HetInfo(uint32_t* ptr) : HetInfo(*ptr, *(ptr+1), *(ptr+2), *(float*)(ptr+3)) {}
	HetInfo() : vcf_line(0), a0(0), a1(0), pp(NAN) {}
	HetInfo(int vcf_line, int a0, int a1, float pp) : vcf_line(vcf_line), a0(a0), a1(a1), pp(pp) {}
	int vcf_line;
	int a0;
	int a1;
	float pp;

	std::string to_string() const {
		std::string output("Position : ");
		output += std::to_string(vcf_line);
		output += std::string(" ") + std::to_string(bcf_gt_allele(a0)) + (bcf_gt_is_phased(a1) ? "|" : "/") + std::to_string(bcf_gt_allele(a1));
		output += std::string(" PP : ") + std::to_string(pp);

		return output;
	}

	void to_stream(std::fstream& ofs) const {
		ofs.write(reinterpret_cast<const char*>(&vcf_line), sizeof(vcf_line));
		ofs.write(reinterpret_cast<const char*>(&a0), sizeof(a0));
		ofs.write(reinterpret_cast<const char*>(&a1), sizeof(a1));
		ofs.write(reinterpret_cast<const char*>(&pp), sizeof(pp));
	}

	void from_stream(std::ifstream& ifs) {
		ifs.read(reinterpret_cast<char*>(&vcf_line), sizeof(vcf_line));
		ifs.read(reinterpret_cast<char*>(&a0), sizeof(a0));
		ifs.read(reinterpret_cast<char*>(&a1), sizeof(a1));
		ifs.read(reinterpret_cast<char*>(&pp), sizeof(pp));
	}
};

constexpr bool operator==(const HetInfo& lhs, const HetInfo& rhs) {
	return lhs.vcf_line == rhs.vcf_line &&
	       lhs.a0 == rhs.a0 && lhs.a1 == rhs.a1 &&
	       ((std::isnan(lhs.pp) && std::isnan(rhs.pp)) || (lhs.pp == rhs.pp));
}

constexpr bool operator!=(const HetInfo& lhs, const HetInfo& rhs) {
	return !(lhs == rhs);
}

class SampleBlock {
public:
	static void write_to_stream(std::fstream& ofs, const std::vector<HetInfo>& his, uint32_t id) {
		const uint32_t mark = 0xd00dc0de;
		const uint32_t size = his.size();
		// Write a mark
		ofs.write(reinterpret_cast<const char*>(&mark), sizeof(uint32_t));
		// Write the ID
		ofs.write(reinterpret_cast<const char*>(&id), sizeof(uint32_t));
		// Write the size
		ofs.write(reinterpret_cast<const char*>(&size), sizeof(uint32_t));
		// Write the Het Infos
		for (auto& hi : his) {
			hi.to_stream(ofs);
		}
	}
};

#endif /* __HET_INFO_HPP__ */

constexpr size_t PLOIDY_2 = 2;

class PPPred {
public:
	PPPred(float pp_threshold) : pp_threshold(pp_threshold) {}

	bool operator()(HetInfo hi) const { return !(std::isnan(hi.pp)) && hi.pp < pp_threshold; }

protected:
	float pp_threshold;
};

template <typename T, class Pred>
class GenericKeepFifo {
public:
	GenericKeepFifo(const size_t size, Pred p) :
	size(size),
	mid(size/2),
	p(p) {
		if (!(size & 1)) {
			std::cerr << "FIFO size should be odd ! Adjusting size to " << ++this->size << std::endl;
			this->mid = this->size/2;
		}
	}

	void insert(T item) {
		// If the FIFO is empty, fill with "dummy items" (to simplify logic)
		if (items.size() == 0) {
			for (size_t i = 0; i < size; ++i) {
				items.push_back(FIFOItem(item));
				// Setting this true will make sure "keep" doesn't save them
				items.back().kept = true;
			}
		}

		// Push item at the end
		items.push_back(FIFOItem(item));

		// Pop the oldest item
		items.pop_front();

		// If predicate (e.g., small PP)
		if (p(items[mid].item)) {
			// Keep the information
			keep();
		}
	}

	void finalize() {
		// Search for predicate (e.g., small PP) at the end
		for (size_t i = mid + 1; i < items.size(); ++i) {
			if (p(items[i].item)) {
				keep_end(i-mid);
				break;
			}
		}
	}

	size_t get_number_kept_with_pred() const {
		size_t n = 0;
		for (auto& i : kept_items) {
			if (p(i)) {
				n++;
			}
		}
		return n;
	}

	const std::vector<T>& get_kept_items_ref() const {
		return kept_items;
	}

private:
	class FIFOItem {
	public:
		FIFOItem(T item) : item(item), kept(false) {}
		T item;
		bool kept;
	};

	inline void keep_end(size_t start) {
		for (size_t i = start; i < items.size(); ++i) {
			if (!items[i].kept) {
				kept_items.push_back(items[i].item);
				items[i].kept = true;
			}
		}
	}

	inline void keep() {
		for (size_t i = 0; i < items.size(); ++i) {
			if (!items[i].kept) {
				kept_items.push_back(items[i].item);
				items[i].kept = true;
			}
		}
	}

protected:
	size_t size;
	size_t mid;
	/// @brief FIFO Items
	std::deque<FIFOItem> items;
	/// @brief Kept items surrounding predicate
	std::vector<T> kept_items;
	/// @brief Predicate that tells us when to keep
	Pred p;
};

class binary2sapphire {
public:

	//PARAM
	std::string region;
	int nthreads;

	//CONSTRUCTORS/DESCTRUCTORS
	binary2sapphire(std::string _region, int _nthreads, float maf_threshold = 0.001,
			bool line_from_vcf = false, size_t fifo_size = 5,
			bool pp_from_maf = false, bool pp_from_af = true);
	~binary2sapphire();

	//PROCESS
	void convert(std::string, std::string);

protected:

	void set_progress(const size_t progress) {
		this->progress = progress;
	}

	void set_maf_threshold(const float maf_threhsold) {
		//std::cout << "Setting MAF threshold to " << maf_threhsold << std::endl;
		MAF_THRESHOLD = maf_threhsold;
	}

	void finalize() {
		// Finalize FIFOs
		for (auto& f : fifos) {
			f.finalize();
		}
	}

	void show_info() {
		size_t total_kept = 0;
		size_t total_kept_pred = 0;
		for (auto& f : fifos) {
			auto kept_items = f.get_kept_items_ref();
			total_kept += kept_items.size();
			total_kept_pred += f.get_number_kept_with_pred();
		}

		vrb.bullet("Extracted a total of " + stb.str(total_kept) + " genotypes");
		vrb.bullet("From which a total of " + stb.str(total_kept_pred) + " were selected given the predicate");
	}

	void write_to_file(std::string filename) {
		std::fstream ofs(filename, std::ios_base::binary | std::ios_base::out | std::ios_base::trunc);
		if (!ofs.is_open()) {
			vrb.error("Cannot open file " + stb.str(filename));
		}

		const uint32_t endianness = 0xaabbccdd;
		// Write endianness
		ofs.write(reinterpret_cast<const char*>(&endianness), sizeof(uint32_t));
		// Write number of samples
		uint32_t num_samples = stop_id-start_id;
		ofs.write(reinterpret_cast<const char*>(&num_samples), sizeof(uint32_t));
		// Write offset table
		uint64_t dummy_offset = 0xdeadc0dedeadc0de;
		auto table_seek = ofs.tellp();
		for (size_t i = start_id; i < stop_id; ++i) {
			ofs.write(reinterpret_cast<const char*>(&dummy_offset), sizeof(uint64_t));
		}

		// Write the data for all the samples
		std::vector<uint64_t> offset_table(stop_id-start_id);
		for (size_t i = start_id; i < stop_id; ++i) {
			const size_t idx = i-start_id;
			offset_table[idx] = ofs.tellp();
			SampleBlock::write_to_stream(ofs, fifos[idx].get_kept_items_ref(), i);
		}

		// Rewrite the offset table
		ofs.seekp(table_seek);
		for (size_t i = start_id; i < stop_id; ++i) {
			const size_t idx = i-start_id;
			ofs.write(reinterpret_cast<const char*>(&offset_table[idx]), sizeof(decltype(offset_table.front())));
		}

		vrb.bullet("Done writing file " + stb.str(filename));

		ofs.close();
	}

	const size_t FIFO_SIZE;
	const float PP_THRESHOLD;
	float MAF_THRESHOLD;
	float *pp_arr;
	int pp_arr_size;
	size_t start_id;
	size_t stop_id;
	size_t line_counter;
	size_t print_counter;
	const PPPred pred;
	size_t progress;
	bool pp_from_maf;
	bool pp_from_af;
	bool extract_acan;
	bool line_from_vcf;
	std::vector<uint32_t> number_of_het_sites;
	std::vector<uint32_t> number_of_low_pp_sites;
	std::vector<uint32_t> number_of_snp_low_pp_sites;
	std::vector<uint32_t> number_of_non_snp;
	std::vector<GenericKeepFifo<HetInfo, PPPred> > fifos;
	std::string search_in_file;
};

#endif
