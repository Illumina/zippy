# ZIPPY
### The
### ZIppy
### Pipeline
### Prototyping
### sYstem

ZIPPY is a powerful, easy-to-use NGS pipeline prototyping system, with batteries included.  ZIPPY helps you create JSON-based pipeline specifications, which it then uses to execute a series of pipeline stages, with no user-written code required.

## With ZIPPY, you can:
- Generate a new pipeline without any coding required
- Auto-generate parameters files from just a list of stages
- Use our ultra-modular system for adding new modules, which can then interface with every other module 
- Make a custom pipeline in 10 minutes or less

## ZIPPY includes a variety of modules, including:
- Bcl2Fastq
- BWA
- Picard alignment stats
- RSEM
- MACS
- BAM subsampling

There will be many more to come!  (See the full list of modules [here](https://github.com/Illumina/zippy/wiki/Zippy-modules)).

## Limitations:
ZIPPY uses black magic that relies upon the CPython implementation of Python.  Currently, only CPython 2.7 is supported.  ZIPPY also requires several python modules.  To make life easier, an executable version of ZIPPY is available (see the releases page!).

### Running ZIPPY from source
If you would like to run ZIPPY from source, there are a couple of things to note.
- You must use CPython 2.7
- You must install the modules 'commentjson' and 'pyflow' (note: pyflow is not in pypi, but can be found [here](https://github.com/Illumina/pyflow)).  You may optionally install the package 'meta' which may improve the coverage of parameters from make_params if you do not have the source of some of your imported modules (e.g., only .pyc files).
- Run the tests file make_params_test.py and see if it's OK!

# Using ZIPPY:
0. Install ZIPPY by using 'pip install zippy-pipeline'

1. Make a proto-workflow file.

    A proto-workflow is a very simple json file that lists the pipeline stages you want, in the order you want them.  Here's a simple proto-workflow that I use for some RNA-seq applications:
    ```
    {
        "stages": ["bcl2fastq", "rsem"]
    }
    ```
Yeah, that's really it.

2. Compile the proto-workflow

    execute 'python -m zippy.make_params my_proto.json my_params.json'

3. Fill in the blanks and/or connect the dots

    Open up my_params.json and fill in all the parameters that are currently blank.

4. Run ZIPPY

    To run ZIPPY, execute 'python -m zippy.zippy my_params.json'


**More information is on the git wiki.**

v2.1.0 (4/11/2018)
- First public release
- Parameters 2.0.  Support for subworkflows, parameter nesting and environments.
- Support for running docker containers through singularity
- Better optional parameters support.  You can now create the function define_optionals(self), which returns a map of default values for optional parameters.
- New stages (Nirvana variant annotation, minimap2 aligner, primer dimer detection)
- New unit tests
- Fewer known bugs

v2.0.0 (2/14/2018)

It's here!  Release 2.0 inaugurates ZIPPY-as-a-package.  You can now pip install ZIPPY, and once you do so, run it as python -m zippy.zippy.  Furthermore, you can import zippy from anywhere and use the API.  Special thanks to Wilfred Li for this work.  As a bonus, we're finally moving to semantic-ish versioning.

- ZIPPY as a package
- API interface should be finalized
- better docker support (can now support docker pull)
- several small bugfixes
- small documentation improvements

v1.9{6} (12/7/2017)
- Provisional ZIPPY API (function names are interfaces not final)
- Removed several external dependencies
- Parameter detection now defaults to the inspect module, instead of the meta module (i.e., we avoid decompilation when possible)
- Support for running locally instead of on SGE
- Memory/core requirements can now be set both per-stage and globally
- New stages (allowing you to delete intermediate outputs, to convert a bam to a fastq, and to downsample in style with a bloom filter)
- DataRunner improvements: can now load sample list from a samplesheet, and can now manually map file names to zippy types
- Better support for modules in external files (make_params now supports them fully)
- Yet fewer bugs!  Or at least... different bugs.

v1.99999 (8/30/17)
- Support for external modules (Runners defined in external code)
- Lots of workflow improvements (argument passing for nearly everything, new bcl2fastq and BWA workflows)
- New workflows (DNA/RNA QC, edger and Salmon)
- New help system (run 'python zippy.py --help')

v1.9999 (5/1/17)
- Unit tests
- Wildcards in the parameters files
- Merge stages
- Support for optional parameters

v1.999 (2/16/17)
- Arbitrary stage chaining
- More stages
- Fewer bugs

v1.99
First end-to-end version of ZIPPY

## License
Copyright 2018 Illumina

Licensed under the Apache License, Version 2.0 (the "License"); you may not use this file except in compliance with the License. You may obtain a copy of the License at http://www.apache.org/licenses/LICENSE-2.0

Unless required by applicable law or agreed to in writing, software distributed under the License is distributed on an "AS IS" BASIS, WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied. See the License for the specific language governing permissions and limitations under the License.
