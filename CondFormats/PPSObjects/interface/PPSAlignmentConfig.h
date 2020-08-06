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
#include <iostream>

//---------------------------------------------------------------------------------------------

struct SelectionRange
{
    double x_min;
    double x_max;

    SelectionRange(double x_min = 0., double x_max = 0.);

    COND_SERIALIZABLE;
};

//---------------------------------------------------------------------------------------------

struct SectorConfig
{
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

//---------------------------------------------------------------------------------------------

class PPSAlignmentConfig
{
public:
    // Getters
    unsigned int fill() const;
    unsigned int xangle() const;
    double beta() const;
    std::string dataset() const;

	std::map<unsigned int, std::string> rpTags() const;

	std::vector<std::string> inputFiles() const;

	std::map<unsigned int, double> alignmentCorrectionsX() const;
    std::map<unsigned int, double> alignmentCorrectionsY() const;

	bool aligned() const;

	double n_si() const;

	SectorConfig sectorConfig45() const;
    SectorConfig sectorConfig56() const;

	std::vector<std::string> matchingReferenceDatasets() const;
	std::map<unsigned int, SelectionRange> matchingShiftRanges() const;

	std::map<unsigned int, SelectionRange> alignment_x_meth_o_ranges() const;
	std::map<unsigned int, SelectionRange> alignment_x_relative_ranges() const;

	std::map<unsigned int, SelectionRange> alignment_y_ranges() const;

    // Setters 
    void setFill(unsigned int fill);
    void setXangle(unsigned int xangle);
    void setBeta(double beta);
    void setDataset(std::string &dataset);

    void setRpTags(std::map<unsigned int, std::string> &rpTags);

    void setInputFiles(std::vector<std::string> &inputFiles);

    void setAlignmentCorrectionsX(std::map<unsigned int, double> &alignmentCorrectionsX);
    void setAlignmentCorrectionsY(std::map<unsigned int, double> &alignmentCorrectionsY);

    void setAligned(bool aligned);

    void setN_si(double n_si);

    void setSectorConfig45(SectorConfig &sectorConfig45);
    void setSectorConfig56(SectorConfig &sectorConfig56);

    void setMatchingReferenceDatasets(std::vector<std::string> &matchingReferenceDatasets);
    void setMatchingShiftRanges(std::map<unsigned int, SelectionRange> &matchingShiftRanges);

    void setAlignment_x_meth_o_ranges(std::map<unsigned int, SelectionRange> &alignment_x_meth_o_ranges);
    void setAlignment_x_relative_ranges(std::map<unsigned int, SelectionRange> &alignment_x_relative_ranges);

    void setAlignment_y_ranges(std::map<unsigned int, SelectionRange> &alignment_y_ranges);
    
    // << operator
    friend std::ostream &operator<<(std::ostream &os, PPSAlignmentConfig c);

private:
    unsigned int fill_;
	unsigned int xangle_;
	double beta_;
	std::string dataset_;

	std::map<unsigned int, std::string> rpTags_;

	std::vector<std::string> inputFiles_;

	std::map<unsigned int, double> alignmentCorrectionsX_, alignmentCorrectionsY_;

	bool aligned_;

	double n_si_;

	SectorConfig sectorConfig45_, sectorConfig56_;

	std::vector<std::string> matchingReferenceDatasets_;
	std::map<unsigned int, SelectionRange> matchingShiftRanges_;

	std::map<unsigned int, SelectionRange> alignment_x_meth_o_ranges_;
	std::map<unsigned int, SelectionRange> alignment_x_relative_ranges_;

	std::map<unsigned int, SelectionRange> alignment_y_ranges_;

    COND_SERIALIZABLE;
};

#endif 