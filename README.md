# bHub\_sim
A repository for simulating islets of Langerhans beta-cells hub phenomenon. This repository supports the results published in <https://doi.org/10.1080/19382014.2018.1493316>.

## Requirements

- Python 2.7 or Python 3.4+
- Python libraries: `numpy` `matplotlib` `neuron`

`neuron` can be downloaded and installed [here](https://www.neuron.yale.edu/neuron/download).

## How to run simulation

This repository is organised this the following way:

- [models](./models) contains all the models used in the study, written in `neuron` readable format. It includes beta cell models and gap junction models.
- [morphologies](./morphologies) contains the x, y, z Cartesian coordinates of the cells in islets of Langerhans, including mouse and human.
- [src](./src), [src-local](./src-local), [src-mpi](./src-mpi), [src-test-pcaer](./src-test-pcaer), [src-delta-cell-sim](./src-delta-cell-sim) contain the 'run' scripts for the simulations.
- [test\_package](./test_package) contains some quick tests (small examples) of the `neuron` package.
- [output](./output) is a placeholder for the outputs (raw traces of the simulated membrane potential and intracellular calcium) from the 'run' scripts in the `src`s.
