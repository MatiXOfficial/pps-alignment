/****************************************************************************
 *
 *  CondFormats/PPSObjects/interface/PPSAlignmentConfig.h
 *
 *  Description : Alignment parameters
 *
 * Authors:
 *  - Jan Kašpar
 *  - Mateusz Kocot
 *
 ****************************************************************************/

#ifndef CondFormats_PPSObjects_PPSAlignmentConfig_h
#define CondFormats_PPSObjects_PPSAlignmentConfig_h

#include "CondFormats/Serialization/interface/Serializable.h"

#include <vector>
#include <string>
#include <map>

//---------------------------------------------------------------------------------------------

struct SelectionRange
{
	double x_min;
	double x_max;

	COND_SERIALIZABLE;
};

//---------------------------------------------------------------------------------------------

struct RPConfig
{
	std::string name;
	unsigned int id;
	std::string position;
	double slope;
	double sh_x;

	double y_max_fit;
	double y_cen_add;
	double y_width_mult;

	COND_SERIALIZABLE;
};
std::ostream &operator<<(std::ostream &os, RPConfig &rc);

//---------------------------------------------------------------------------------------------

struct SectorConfig
{
	std::string name;
	RPConfig rp_N, rp_F;
	double slope;

	bool cut_h_apply;
	double cut_h_a, cut_h_c, cut_h_si;

	bool cut_v_apply;
	double cut_v_a, cut_v_c, cut_v_si;

	int nr_x_slice_n;
	double nr_x_slice_min, nr_x_slice_w;
	int fr_x_slice_n;
	double fr_x_slice_min, fr_x_slice_w;

	COND_SERIALIZABLE;
};
std::ostream &operator<<(std::ostream &os, SectorConfig &sc);

//---------------------------------------------------------------------------------------------

struct Binning
{
	double bin_size_x;	// mm
	unsigned int n_bins_x;

	double pixel_x_offset;

	unsigned int n_bins_y;
	double y_min, y_max;
};
std::ostream &operator<<(std::ostream &os, Binning &b);

//---------------------------------------------------------------------------------------------

class PPSAlignmentConfig
{
public:
	// Getters
	std::vector<std::string> sequence() const;

	SectorConfig sectorConfig45() const;
	SectorConfig sectorConfig56() const;

	std::map<unsigned int, double> alignmentCorrectionsX() const;
	std::map<unsigned int, double> alignmentCorrectionsY() const;

	bool aligned() const;
	double x_ali_sh_step() const;

	double y_mode_sys_unc() const;
	double chiSqThreshold() const;
	double y_mode_unc_max_valid() const;
	double y_mode_max_valid() const;

	double n_si() const;

	std::vector<std::string> matchingReferenceDatasets() const;
	std::map<unsigned int, SelectionRange> matchingShiftRanges() const;

	std::map<unsigned int, SelectionRange> alignment_x_meth_o_ranges() const;
	std::map<unsigned int, SelectionRange> alignment_x_relative_ranges() const;

	std::map<unsigned int, SelectionRange> alignment_y_ranges() const;

	Binning binning() const;

	// Setters 
	void setSequence(std::vector<std::string> &sequence);

	void setSectorConfig45(SectorConfig &sectorConfig45);
	void setSectorConfig56(SectorConfig &sectorConfig56);

	void setAlignmentCorrectionsX(std::map<unsigned int, double> &alignmentCorrectionsX);
	void setAlignmentCorrectionsY(std::map<unsigned int, double> &alignmentCorrectionsY);

	void setAligned(bool aligned);
	void setX_ali_sh_step(double x_ali_sh_step);

	void setY_mode_sys_unc(double y_mode_sys_unc);
	void setChiSqThreshold(double chiSqThreshold);
	void setY_mode_unc_max_valid(double y_mode_unc_max_valid);
	void setY_mode_max_valid(double y_mode_max_valid);

	void setN_si(double n_si);

	void setMatchingReferenceDatasets(std::vector<std::string> &matchingReferenceDatasets);
	void setMatchingShiftRanges(std::map<unsigned int, SelectionRange> &matchingShiftRanges);

	void setAlignment_x_meth_o_ranges(std::map<unsigned int, SelectionRange> &alignment_x_meth_o_ranges);
	void setAlignment_x_relative_ranges(std::map<unsigned int, SelectionRange> &alignment_x_relative_ranges);

	void setAlignment_y_ranges(std::map<unsigned int, SelectionRange> &alignment_y_ranges);

	void setBinning(Binning &binning);

	// << operator
	friend std::ostream &operator<<(std::ostream &os, PPSAlignmentConfig c);

private:
	std::vector<std::string> sequence_;

	SectorConfig sectorConfig45_, sectorConfig56_;

	std::map<unsigned int, double> alignmentCorrectionsX_, alignmentCorrectionsY_;

	bool aligned_;
	double x_ali_sh_step_;	// mm

	double y_mode_sys_unc_;
	double chiSqThreshold_;
	double y_mode_unc_max_valid_;
	double y_mode_max_valid_;

	double n_si_;

	std::vector<std::string> matchingReferenceDatasets_;
	std::map<unsigned int, SelectionRange> matchingShiftRanges_;

	std::map<unsigned int, SelectionRange> alignment_x_meth_o_ranges_;
	std::map<unsigned int, SelectionRange> alignment_x_relative_ranges_;

	std::map<unsigned int, SelectionRange> alignment_y_ranges_;

	Binning binning_;

	COND_SERIALIZABLE;
};

std::ostream &operator<<(std::ostream &os, PPSAlignmentConfig c);

#endif 