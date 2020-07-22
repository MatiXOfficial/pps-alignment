# pps-alignment

## Installation
```
scram project CMSSW_11_2_X_2020-07-16-1100
cd CMSSW_11_2_X_2020-07-16-1100/src
cmsenv
git clone https://github.com/MatiXOfficial/pps-alignment.git .
scram b
```
## Event Setup example (does not work)
```
cd CalibPPS/ESProducers/test
cmsRun test_PPSAlignmentConfigESSource.py
```
## Relevant files
- DQM plugins
  - `CalibPPS/Alignment/plugins`
- Event setup class
  - `CondFormats/PPSObjects/interface/PPSAlignmentConfig.h`
- ESProducer
  - `CalibPPS/ESProducers/plugins/PPSAlignmentConfigESSource.cc`
- ES Record
  - `CondFormats/DataRecord/interface/PPSAlignmentConfigRcd.h`
  - `CondFormats/DataRecord/src/PPSAlignmentConfigRcd.cc`