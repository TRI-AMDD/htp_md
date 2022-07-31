# htp_md
htp_md is the analysis module in a suite of tools ([UI](htpmd.matr.io)) that streamlines the process of analyzing, storing, visualizing, and predicting properties based on raw trajectory data, in support of research by polymers program at Massachusetts Institute of Technology and the Toyota Research Institute.

htp_md extracts the following properties from raw trajectory data:
* SMILES (input metadata)
* Simulation length (input metadata)
* Ionic Conductivity
* Diffusion Coefficient (ions and polymer chains)
* Molality
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

4. Start the container and the environment. If you would like to use provided test data, use the following command:
```
docker run -it htpmd
conda activate htpmd
```
If you would instead like to use your own data, use the following command:
```
docker run -it -v full/path/to/your/own/data:/src/your_data htpmd
conda activate htpmd
```
where `full/path/to/your/own/data` is on your local computer, and `your_data` is the destination folder within the container.

Once in the container, use commands below to run analysis functions.

### Test Data
We package some data which can be used for testing purposes. This data can be found in `./test_data/`, including several datasets:
- `nacl_water`: a trajectory of a 1m aqueous NaCl electrolyte at 350 K. The system is composed of 222 water molecules and 20 ion pairs. The trajectory is 2 ns long, with 201 snapshots. The interatomic potential used is the SPC/E model for water, with standard Joung-Cheatham parameters for the ions. All parameters are reported in [`10.1021/jp902584c`](https://doi.org/10.1021/jp902584c).
- `9-0-246295613-0`: a small fraction from a trajectory of polymer electrolytes with LiTFSI at 353 K. The trajectory is 14 ps long, with 7 snapsshots. The interatomic potential used is PCFF+, with the charge distribution of TFSI- adjusted. Details for the simualtion can be found in [`10.1038/s41467-022-30994-1`](https://doi.org/10.1038/s41467-022-30994-1).
- `9-0-413610210-0`: a small fraction from a trajectory of polymer electrolytes with LiTFSI at 353 K. The trajectory is 14 ps long, with 8 snapsshots. The interatomic potential used is PCFF+, with the charge distribution of TFSI- adjusted. Details for the simualtion can be found in [`10.1038/s41467-022-30994-1`](https://doi.org/10.1038/s41467-022-30994-1).

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

To run analysis functions on your own volume attached to the container (following instructions above):
```
python main.py <action> [-d </src/your/data>]
```

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
Contact the HTP team at Toyota Research Institute (materials-support@tri.global) with your name, affiliation, and a description of your data. Any contributed data should be reproducible and are required to include the following:
* `Metadata.json`, which contains the description of the system, such as SMILES string for monomer and polymer, force field used, material group (polymer), temperature, time step, and cation and anion information (name and atom type), as well as polymer information (range of atom types)
* `relaxed.lmp`, which contains the starting configuration for LAMMPS production run, in the LAMMPS data file format (https://docs.lammps.org/2001/data_format.html)
* `in.lmp`, which contains input information for running LAMMPS, in the LAMMPS input script format (https://docs.lammps.org/Commands_input.html)

## Authors
Toyota Research Institute
- Ha-Kung Kwon
- Daniel Schweigert
- Arash Khajeh

Massachusetts Institute of Technology
- Tian Xie
- Arthur France-Lanord
- Emily Crabb
- Sheng Gong

## How to Cite
If you use htp_md, please cite the following: TODO
