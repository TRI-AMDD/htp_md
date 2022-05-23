# htp_md
htp_md is the analysis module in a suite of tools ([htp_md_worker],[UI]) that streamlines the process of analyzing, storing, visualizing, and predicting properties based on raw trajectory data, in support of research by polymers program at Massachusetts Institute of Technology and the Toyota Research Institute.

Facilitated by [htp_md_worker](https://github.com/tri-amdd/htp_md_worker), htp_md extracts the following properties from raw trajectory data:
* Molality (input metadata)
* SMILES (input metadata)
* Simulation length (input metadata) 
* Li-ion Conductivity
* Diffusion Coefficient (ions and polymer chains)
* Transference Number
* Structure (CIF)
* Mean Squared Displacement (MSD) time series of ions

## Getting Started
### Dependencies
Dependencies are found in `requirements.txt`. 

### Installation
1. Clone the repo and install.
```
git clone git@github.com:TRI-AMDD/htp_md.git
```
2. Build the docker container
```
docker build -t htpmd .
```
3. To test that the build is complete and runs properly, run pytests:
```
docker run htpmd python -m pytest
```

4. Start the container and the environment with the following command:
```
docker run -it htpmd
conda activate htpmd
```
Once in the container, use commands below to run analysis functions. 

### Test Data 
We package some data which can be used for testing purposes. This data can be found in `./test_data/`, including several datasets: 
- `nacl_water`: a trajectory of a 1m aqueous NaCl electrolyte at 350 K. The system is composed of 222 water molecules and 20 ion pairs. The trajectory is 2 ns long, with 201 snapshots. The interatomic potential used is the SPC/E model for water, with standard Joung-Cheatham parameters for the ions. All parameters are reported in `10.1021/jp902584c`. 
- `9-0-246295613-0`: a small fraction from a trajectory of polymer electrolytes with LiTFSI at 353 K. The trajectory is 14 ps long, with 7 snapsshots. The interatomic potential used is PCFF+, with the charge distribution of TFSI- adjusted. Details for the simualtion can be found in arXiv:2101.05339. 
- `9-0-413610210-0`: a small fraction from a trajectory of polymer electrolytes with LiTFSI at 353 K. The trajectory is 14 ps long, with 8 snapsshots. The interatomic potential used is PCFF+, with the charge distribution of TFSI- adjusted f
ollowing 10.1021/jp077026y. Details for the simualtion can be found in arXiv:2101.05339.

## Using the functions
To get a dictionary of all results:

```python
import htpmd
results = htpmd.analyze('test_data/9-0-246295613-0')
```
This returns a dictionary of results. To see subfields of results
`results.keys()`


### Testing
To run unit tests, run:
`python -m pytest` 

## Using htp_md
To run analysis functions, run:
`python main.py <action> [-d <dir_path>]`

## How to contribute
User contributions for new analysis functions and data are greatly appreciated. 

### Contributing a new analysis function
1. Fork the Project
2. Create your Feature Branch (`git checkout -b feature/AmazingFeature`)
3. Commit your Changes (`git commit -m 'Add AmazingFeature'`)
4. Push to the Branch (`git push origin feature/AmazingFeature`)
5. Open a Pull Request

When contributing a new function, please follow the [template](https://github.com/TRI-AMDD/htp_md/blob/master/src/htpmd/shared/template.py). Each Pull Request for a new function should contain the following:
* `function.py`, following template guidelines
* `function_test.py`, following template guidelines
* Test data and expected outcome

### Contributing a new trajectory
Contact the HTP team at Toyota Research Institute (materials-support@tri.global) with your name, affiliation, and a description of your data. 

## Version History
* 0.2.4
    * Bug fixes

## Authors
Toyota Research Institute
Massachusetts Institute of Technology (Tian, Sheng, Arthur, Emily)

## How to Cite
If you use htp_md, please cite the following: TODO
