# pps-alignment

## Installation
```
scram project CMSSW_11_2_0_pre2
cd CMSSW_11_2_0_pre2/src
cmsenv
git clone https://github.com/MatiXOfficial/pps-alignment.git .
scram b -j 8
```
## Event Setup example
```
cd CalibPPS/ESProducers/test
cmsRun test_PPSAlignmentConfigESSource_cfg.py
```
<!-- ## Relevant files
- DQM plugins
  - `CalibPPS/Alignment/plugins`
- Event setup class
  - `CondFormats/PPSObjects/interface/PPSAlignmentConfig.h`
- ESProducer
  - `CalibPPS/ESProducers/plugins/PPSAlignmentConfigESSource.cc`
- ES Record
  - `CondFormats/DataRecord/interface/PPSAlignmentConfigRcd.h`
  - `CondFormats/DataRecord/src/PPSAlignmentConfigRcd.cc` -->