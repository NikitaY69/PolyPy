# Static Monte Carlo Lattice Polymer simulation program

This package provides a Monte Carlo computational approach to simulate 3-dimensional polymers on cubic lattices.

## Description
Our program runs Monte Carlo simulations to compute 3-dimensional lattice-polymers observables such as length, gyration length and extension. It is based on a static sampling: the **Rosenbluth method**. The sampling can be enhanced thanks to the PERM algorithm. The evaluation of the *error* is based on bootstrapping. Visited-sites heatmaps can help picture high-weight configurations on top of single random walks.  

### What can our program simulate ?
- Self-Avoiding Random Walks
- Interacting Self-Avoiding Walks (nearest non-paired neighbors attraction)
- Biased Self-Avoiding Walks (constant force constraint)

Please refer to the `report/` folder for a brief review of the project's scientific background and our results.

### Example
An example of usage is provided in the `results.py` file. A full how-to tutorial will soon be published as a notebook in the present repository. 

## Installation

### Cloning 
In your favourite terminal, run:
```
git clone https://github.com/NikitaY69/PolyPy.git
cd PolyPy
```

### Requirements
Our project requiring very few specific packages, you can just run 
```
pip install -r requirements.txt
```
to make sure you possess all the prerequisites. In case you want to work on a seperate environment, proceed as follows:
- Create new virtual environment
```
python -m venv polymers
```
- Activate it
```
source polymers/bin/activate
```
- Upgrade pip
```
pip install --upgrade pip
```
- Install required packages
```
pip install -r requirements.txt
```

### Installing editable version
Because our program is not yet in a perfect stance, we recommend to install it with the edit flag in case you want to adapt the code to your specific needs. To do so, you must run the following command in the present directory:
```
pip install -e .
```