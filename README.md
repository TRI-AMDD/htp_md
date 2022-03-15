# htp_md
htp_md is the analysis module in a suite of tools ([htp_md_worker],[UI]) that streamlines the process of analyzing, storing, visualizing, and predicting properties based on raw trajectory data, in support of research by polymers program at Massachusetts Institute of Technology and the Toyota Research Institute.

Facilitated by [htp_md_worker](https://github.com/tri-amdd/htp_md_worker), htp_md extracts the following properties from raw trajectory data:
* Molality (input metadata)
* SMILES (input metadata)
* Conductivity
* Lithium diffusivity
* Anion diffusivity
* Polymer diffusivity
* Transference Number
* Structure (CIF)
* Lithium MSD
* Anion MSD

## Getting Started
### Dependencies
Dependencies are found in `env.yml` and `requirements.txt`. 

### Installation
1. Clone the repo and install. We recommend having the latest version of conda.
```git clone git@github.com:TRI-AMDD/htp_md.git
cd htp_md
conda update conda
conda env create -f env.yml
```

Alternatively, one can create an environment separately and install the packages manually. 

```
conda create -n htpmd python=3.6
conda install --name htpmd --file spec-file.txt
```
Once this environment is installed, activate and continue through set-up.

```
conda activate htpmd
pip install .

```

### Test Data 
We package some data which can be used for testing purposes. This data can be found in `./test_data/`, including several datasets: 
- `nacl_water`: a trajectory of a 1m aqueous NaCl electrolyte at 350 K. The system is composed of 222 water molecules and 20 ion pairs. The trajectory is 2 ns long, with 201 snapshots. The interatomic potential used is the SPC/E model for water, with standard Joung-Cheatham parameters for the ions. All parameters are reported in `10.1021/jp902584c`. 
- `9-0-246295613-0`: ?
- `9-0-413610210-0`: ?

## Using the functions
To get a dictionary of all results:

```python
import htpmd
results = htpmd.analyze('test_data/9-0-246295613-0')
```
This returns a dictionary of results. To see subfields of results
`results.keys()`

TODO: A tutorial that walks through installation of repo, download a test trajectory, use analyze_all to extract properties, use ML to make a prediction, and visualize. 

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

When contributing a new function, please follow the [template](https://github.com/TRI-AMDD/htp_md/blob/master/htpmd/shared/template.py). Each Pull Request for a new function should contain the following:
* `function.py`, following template guidelines
* `function_test.py`, following template guidelines
* Test data and expected outcome

### Contributing a new trajectory
Contact the HTP team at Toyota Research Institute (materials-support@tri.global) with your name, affiliation, and a description of your data. 

## Version History
* 0.1.5
    * First release to the public

## Authors
Toyota Research Institute
Massachusetts Institute of Technology (Tian, Sheng, Arthur, Emily)

## How to Cite
If you use htp_md, please cite the following: TODO
