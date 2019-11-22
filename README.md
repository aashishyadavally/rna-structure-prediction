# RNA Structure Prediction
## Synopsis
This mini-project was implemented as a part of UGA's CSCI 6470 Algorithms course.

The aim of this project is to compute the maximum number of permissible letter pairs from the input sequence of 4 letters (A, C, G, U). For these letters, legal
pairs are A-U, U-A, G-C, C-G.  This problem is called *Structure Prediction* because it is a meaningful abstraction of the significant problem RNA secondary structure
prediction in computational biology.  The letters in the input sequence are replaced with denotational symbols '{' and '}' underlining the sequence indicate
which two letters are predicted to be paired and the dot symbol ’.’ indicates a letter predicted to be not paired with any other.

For example,<br/>
**Input**: AUGUCAGCGUU <br/>
**Output**: { } { . } { { } . } . <br/>

### Reference
Sean Eddy, How do RNA Folding Algorithms Work? Nature Biotechnology,
Vol 22, No. 7, July 2004.

## Getting Started
This section describes the preqrequisites, and contains instructions, to get the project up and running.

### Usage
 The user can get a description of the options by using the command: `$ python __main__.py --help`.

 User Options:<br />
 -  **--energy**        <p>If True, returns specific gamma-scores for each base-
                        pair. If False, returns 1 for each base-pair.</p>
 - **--gap**             <p>Number of bases between two nucleotide bases in
                        sequence.</p>
 - **--include_AU**      <p>If True, includes A-U base pair in DP-Table
                        formulation.</p> 

### Built With
* [Python](https://www.python.org/)

## Contributing Guidelines
There are no specific guidelines for contributing, apart from a few general rules such as:
* Code should follow PEP8 standards as closely as possible
* We use [Google-Style docstrings](https://sphinxcontrib-napoleon.readthedocs.io/en/latest/example_google.html) to document the Python modules in this project.

If you see something that could be improved, send a pull request! 
I am always happy to look at improvements, to ensure that `rna-structure-prediction`, as a project, is the best version of itself. 

If you think something should be done differently (or is just-plain-broken), please create an issue.

## License
See the [LICENSE](https://github.com/aashishyadavally/rna-structure-prediction/blob/master/LICENSE) file for more details.
