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
TODO

### Installation
1. Clone the repo and install.
```git clone git@github.com:TRI-AMDD/htp_md.git
cd htp_md
pip install .
```
### Testing
To run unit tests, run:
`python -m pytest` 

## Using htp_md
To run analysis functions, run:
`python main.py <action> [-d <dir_path>]`

### Tutorial
TODO: A tutorial (doc or jupyter notebook) that walks through installation of repo, download a test trajectory, use analyze_all to extract properties, use ML to make a prediction, and visualize. 

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
Contact Toyota Research Institute with your name, affiliation, and a description of your data. 

## Version History
* 0.1.5
    * First release to the public

## Authors
Toyota Research Institute
Massachusetts Institute of Technology (Tian, Sheng, Arthur, Emily)

## How to Cite
If you use htp_md, please cite the following: TODO