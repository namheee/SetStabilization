# SetStabilization

## Usage

+ You can simply download SetStabilization from this git repository, while setup.py is not provided.
+ We provide a main simulation example using Jupyter notebooks. Please refer to `Tutorial.ipynb`.
+ An example network is provided in the `example` directory.  


## Notes
+ To generate random Boolean networks, we used [BNGenerator](https://github.com/choonlog/OutputStabilization), a software that generates a random Boolean network using Biological Boolean logics extracted from 78 Biological Boolean networks in the Cell Collective (https://cellcollective.org/)
+ To identify all attractors of Boolean networks, we utilized [BooleanSim](https://github.com/jehoons/BooleanSim), a Python 3 tool developed based upon the [booleannet](https://github.com/ialbert/booleannet).
+ Main functions are included in `setMain.py` and `setFunction.py`.
+ An auxiliary function is included in `choonFunction.py`, with a referecne to its source at [AttractorTransition](https://github.com/choonlog/AttractorTransition).


## Reference paper
+ A related paper is expected to be published soon.
