
**WARNING! starting from verion 3.2.0 (Aug 22, 2022), [DENOPTIM](https://github.com/denoptim-project/DENOPTIM) provides all the functionality of GM3DFragmenter in a more user-friendly framework that includes a graphical user interface. A number of other improvements to the fragmentation algorithm are available in the DENOPTIM-integrated implementation. Therefore, GM3DFragmenter is discontinued.**


# GM3DFragmenter

GM3DFragmenter is a molecular fragmentation tool meant to support peculiar chemical features, such as multihapto bonds and situations beyond standard valence rules.

This fragmenter takes chemical objects containing formally defined bonds (i.e., a connectivity table) and generates 2D and 3D fragments according to the principles of the class-based approach described in these papers:

1) Foscato, M.; Occhipinti, G.; Venkatraman, V.; Alsberg, B. K.; Jensen, V. R.; Automated Design of Realistic Organometallic, Molecules from Fragments; *J. Chem. Inf. Model.* **2014**, 54, 767–780.
2) Foscato, M.; Venkatraman, V.; Occhipinti, G.; Alsberg, B. K.; Jensen, V. R.; Automated Building of Organometallic Complexes from 3D Fragments; *J. Chem. Inf. Model.* **2014**, 54, 1919–1931.

## Overview of the folder tree
* <code>data</code>: contains various versions of cutting rules
* <code>doc</code>: contains the user manual and, after running the createDoc.sh, the API documentation
* <code>jar</code>: after building (see below), contains the JAR archive files
* <code>lib</code>: contains the third parties library used by GM3DFragmenter
* <code>src</code>: source code of GM3DFragmenter
* <code>test</code>: input and data required to run the functionality tests
* <code>utils</code>: utilities associated with the generation and management of fragments libraries


## Compile GM3DFragmenter

Building scripts are made available for Unix systems. These scripts assume <code>java</code> and <code>javac</code> commands are in your <code>$PATH</code>. To verify this condition, please see if you get a reasonable output from the following commands.

        java -version
        javac -version

To build GM3DFragmenter run:

        ./build_GM3DFragmenter.sh

Now, you are ready to go. Happy fragmenting!

## User Manual
The user manual is available under the <code>doc</code> folder. An <a href="http://htmlpreview.github.io/?https://github.com/denoptim-project/GM3DFragmenter/blob/master/doc/user_manual.html">online version</a> is also available.

**WARNING! starting from verion 3.2.0 (Aug 22, 2022), [DENOPTIM](https://github.com/denoptim-project/DENOPTIM) provides all the functionality of GM3DFragmenter in a more user-friendly framework that includes a graphical user interface. A number of other improvements to the fragmentation algorithm are available in the DENOPTIM-integrated implementation. Therefore, GM3DFragmenter is discontinued.**

## License
The software is distributed under the GNU Affero General Public License version 3. Modifications of the CDK library (affecting org/openscience/cdk/smiles/smarts/SMARTSQueryTool, org/openscience/cdk/smiles/FixBondOrdersTool, org/openscience/cdk/tools/periodictable/PeriodicTable) are discributed under the GNU Lesser General Public License version 2.1. See under the <code>lib</code>.



