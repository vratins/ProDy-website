.. _prody-basics:

ProDy Basics
===============================================================================


We start with importing everything from ProDy package:

.. ipython:: python

   from prody import *
   from pylab import *
   ion()


Functions and classes are named such that they should not create a conflict
with any other package.  In this part we will familiarize with different
categories of functions and methods.


File Parsers
-------------------------------------------------------------------------------

Let's start with parsing a protein structure and then keep working on that
in this part.  File parser function names are prefixed with ``parse``.
You can get a list of parser functions by pressing :kbd:`TAB` after typing
in ``parse``:

.. ipython::

   @verbatim
   In [2]: parse<TAB>
      parseArray        parseCIFStream    parseEMD          parseHiC          parseMSA          
      parsePDBHeader    parsePQR          parseSTAR         parseChainsList   parseDCD          
      parseEMDStream    parseHiCStream    parseNMD          parsePDBStream    parsePSF          
      parseSTRIDE       parseCIF          parseDSSP         parseHeatmap      parseModes        
      parsePDB          parseMMTF         parsePfamPDBs     parseSparseMatrix  

When using :func:`.parsePDB`, usually an identifier will be sufficient,
If corresponding file is found in the current working directory, it will be
used, otherwise it will be downloaded from PDB servers.

Let's parse structure :pdb:`1p38` of p38 MAP kinase (MAPK):

.. ipython:: python

   p38 = parsePDB('1p38')  # returns an AtomGroup object
   p38 # typing in variable name will give some information

We see that this structure contains 2962 atoms.

Now, similar to listing parser function names, we can use tab completion to
inspect the ``p38`` object:

.. ipython::

   @verbatim
   In [2]: p38.num<TAB>
   p38.numAtoms      p38.numChains     p38.numFragments  p38.numSegments
   p38.numBonds      p38.numCoordsets  p38.numResidues


This action printed a list of methods with `num` prefix. Let's use some of
them to get information on the structure:


.. ipython::

   In [2]: p38.numAtoms()
   Out[2]: 2962

   In [2]: p38.numCoordsets()  # returns number of models
   Out[2]: 1

   In [2]: p38.numResidues()  # water molecules also count as residues
   Out[2]: 480


Analysis Functions
-------------------------------------------------------------------------------

Similar to parsers, analysis function names start with ``calc``:

.. ipython::

   @verbatim
   In [2]: calc<TAB>
      calcADPAxes                   calcChainsNormDistFluct       calcCrossProjection           
      calcDistFlucts                calcFractVariance             calcADPs                      
      calcCollectivity              calcCumulOverlap              calcENM                       
      calcGNM                       calcAngle                     calcCovariance                
      calcDeformVector              calcEnsembleENMs              calcGyradius                     
      calcANM                       calcCovOverlap                calcDihedral                  
      calcEnsembleSpectralOverlaps  calcMeff                      calcCenter                    
      calcCrossCorr                 calcDistance                  calcEntropyTransfer           
      calcMSAOccupancy              calcMSF                       calcPairDeformationDist       
      calcPsi                       calcSignatureCollectivity     calcSpecDimension                
      calcOccupancies               calcPercentIdentities         calcRankorder                 
      calcSignatureCrossCorr        calcSpectralOverlap           calcOmega                     
      calcPerturbResponse           calcRMSD                      calcSignatureFractVariance    
      calcSqFlucts                  calcOverallNetEntropyTransfer calcPhi                       
      calcRMSF                      calcSignatureOverlaps         calcSubspaceOverlap              
      calcOverlap                   calcProjection                calcShannonEntropy            
      calcSignatureSqFlucts         calcTempFactors               calcTransformation
      calcTree                  

Let's read documentation of :func:`.calcGyradius` function and use it to
calculate the radius of gyration of p38 MAPK structure:

.. ipython::

   ? calcGyradius

   calcGyradius(p38)


Plotting Functions
-------------------------------------------------------------------------------

Likewise, plotting function names have ``show`` prefix and here is a list of them:

.. ipython::

   @verbatim
   In [2]: show<TAB>
      showAlignment             showCrossProjection       showDomainBar             
      showHeatmap               showMeanMechStiff         showNormDistFunct          
      showAtomicLines           showCumulFractVars        showDomains               
      showLines                 showMechStiff             showNormedSqFlucts         
      showAtomicMatrix          showCumulOverlap          showEllipsoid             
      showLinkage               showMode                  showOccupancies           
      showContactMap            showDiffMatrix            showEmbedding             
      showMap                   showMSAOccupancy          showOverlap                
      showCrossCorr             showDirectInfoMatrix      showFractVars             
      showMatrix                showMutinfoMatrix         showOverlaps 
      showOverlapTable          showScaledSqFlucts        showSignatureCollectivity 
      showSignatureSqFlucts     showVarianceBar           showPairDeformationDist   
      showSCAMatrix             showSignatureCrossCorr    showSignatureVariances                               
      showPerturbResponse       showShannonEntropy        showSignatureDistribution 
      showSqFlucts              showProjection            showSignature1D           
      showSignatureMode         showTree                  showProtein               
      showSignatureAtomicLines  showSignatureOverlaps     showTree_networkx                                                            


We can use :func:`.showProtein` function to make a quick plot of p38 structure:

.. ipython:: python

   @savefig prody_tutorial_basics_protein.png width=4in
   showProtein(p38);

This of course does not compare to any visualization software that you
might be familiar with, but it comes handy to see what you are dealing with.



Protein Structures
-------------------------------------------------------------------------------

Protein structures (:file:`.pdb` or :file:`.cif` files) will be the standard input for most
*ProDy* calculations, so it is good to familiarize with ways to access and
manage PDB file resources.

Fetching PDB files
^^^^^^^^^^^^^^^^^^

First of all, *ProDy* downloads PDB files when needed (these are compressed on the PDB webserver). 
If you prefer saving decompressed files, you can use :func:`.fetchPDB` function as
follows:

.. ipython:: python

  fetchPDB('1p38', compressed=False)

Note that ProDy functions that fetch files or output files return filename
upon successful completion of the task.  You can use this behavior to
shorten the code you need to write, e.g.:

.. ipython:: python

  parsePDB(fetchPDB('1p38', compressed=False)) # same as p38 parsed above

We downloaded and save an uncompressed PDB file, and parsed it immediately.

PDB file resources
^^^^^^^^^^^^^^^^^^

Secondly, ProDy can manage local mirrors of the PDB server or a local PDB folder,
as well as using a server close to your physical location for downloads:

  * One of the `wwPDB`_ FTP servers in US, Europe or Japan can be picked for
    downloads using :func:`.wwPDBServer`.

  * A local PDB mirror can be set for faster access to files using
    :func:`.pathPDBMirror`.


  * A local folder can be set for storing downloaded files for future access
    using :func:`.pathPDBFolder`.

If you are in the Americas now, you can choose the PDB server in the US
as follows:

.. ipython:: python

   wwPDBServer('us')

If you would like to have a central folder, such as :file:`~/Downloads/pdb`,
for storing downloaded PDB files (you will need to make it), do as follows:

.. ipython::
   :verbatim:

   In [19]: mkdir ~/Downloads/pdb;

   In [20]: pathPDBFolder('~/Downloads/pdb')

Note that when these functions are used, ProDy will save your settings
in :file:`.prodyrc` file stored in your home folder.

.. _wwPDB: http://www.wwpdb.org/

..
  :func:`.parsePDB` function is very flexible and can be extremely
  efficient depending on what you want to extract from a PDB file.  It can be
  used to parse specific chains, models, alternate locations, or well-defined
  subsets of atoms from a file.  A detailed usage example can be found in
  :ref:`parsepdb`.

  ProDy can parse other file types, including :file:`.psf` and :file:`.pqr` files.
  All of the functions for accessing and handling protein structural data are
  described in :mod:`.proteins` module reference documentation.
  Also, :ref:`fetchpdb` and :ref:`blastpdb` examples show other ways to
  access the Protein Data Bank (|pdb|) content.

Atom Groups
-------------------------------------------------------------------------------

As you might have noticed, :func:`.parsePDB` function returns structure data
as an :class:`.AtomGroup` object. Let's see for ``p38`` variable from above:

.. ipython:: python

   p38

You can also parse a list of :file:`.pdb` files into a list of :class:`.AtomGroup`
objects:

.. ipython:: python

   ags = parsePDB('1p38', '3h5v')
   ags

If you want to provide a list object you need to provide an asterisk (``*``) to 
let Python know this is a set of input arguments:

.. ipython:: python

   pdb_ids = ['1p38', '3h5v']
   ags = parsePDB(pdb_ids)
   ags 


Data from this object can be retrieved using ``get`` methods. For example:

.. ipython:: python

   p38.getResnames()
   p38.getCoords()


To get a list of all methods use tab completion, i.e. ``p38.<TAB>``.
We will learn more about atom groups in the following chapters.

Indexing
^^^^^^^^

An individual :class:`.Atom` can be accessed by indexing :class:`.AtomGroup`
objects:

.. ipython:: python

   atom = p38[0]
   atom

Note that all ``get/set`` functions defined for :class:`.AtomGroup`
instances are also defined for :class:`.Atom` instances, using singular
form of the function name.

.. ipython:: python

   atom.getResname()

Slicing
^^^^^^^

It is also possible to get a slice of an :class:`.AtomGroup`. For example,
we can get every other atom as follows:

.. ipython:: python

   p38[::2]

Or, we can get the first 10 atoms, as follows:


.. ipython:: python

   p38[:10]

Hierarchical view
^^^^^^^^^^^^^^^^^

You can also access specific chains or residues in an atom group.  Indexing
by a single letter identifier will return a :class:`.Chain` instance:

.. ipython:: python

   p38['A']

Indexing atom group with a chain identifier and a residue number will return
:class:`.Residue` instance:


.. ipython:: python

   p38['A', 100]


See :ref:`atomic` for details of indexing atom groups and :ref:`hierview`
for more on hierarchical views.



ProDy Verbosity
-------------------------------------------------------------------------------

Finally, you might have noticed that ProDy prints some information to the
console after parsing a file or doing some calculations. For example, PDB
parser will print what was parsed and how long it took to the screen::

  @> 1p38 (./1p38.pdb.gz) is found in the target directory.
  @> 2962 atoms and 1 coordinate sets were parsed in 0.08s.

This behavior is useful in interactive sessions, but may be problematic for
automated tasks as the messages are printed to stderr.  The level of verbosity
can be controlled using :func:`.confProDy` function, and calling it as
``confProDy(verbosity='none')`` will stop all information messages permanently.

Reading in Hadoop Sequence Files using PySpark
-------------------------------------------------------------------------------

The Macromolecular Transmission Format (MMTF) is a compact binary format to transmit and store biomolecular structures (https://mmtf.rcsb.org). MMTF files or structure objects can be read into ProDy as an atom group by:

.. ipython:: python

  parseMMTF('1p38')
  parseMMTF('1p38.mmtf')
  parseMMTF(1p38_mmtf_structure_object)

For large scale calculations with distributed parallel frameworks such as PySpark we recommend the use of Hadoop Sequence Files. The entire PDB archive can be downloaded in the full and reduced representations:
https://mmtf.rcsb.org/v1.0/hadoopfiles/full.tar and https://mmtf.rcsb.org/v1.0/hadoopfiles/reduced.tar


We can download the hadoop sequence files using wget or curl:

.. code:: ipython3

    $ wget https://mmtf.rcsb.org/v1.0/hadoopfiles/full.tar
    $ tar -xvf full.tar

or

.. code:: ipython3

    $ curl -O https://mmtf.rcsb.org/v1.0/hadoopfiles/full.tar
    $ tar -xvf full.tar

This will download and unpack the content of the Hadoop Sequence File to
a folder.

In order to perform any analysis using PySpark
(https://spark.apache.org/docs/latest/api/python/index.html), we will
need initizalize a Spark session. Make sure that the proper packages are
imported.

You can initialize PySpark with certain parameters to optimize
perfomance, and adjust based on the available specifications of your
computing environments.

**Driver memory**: This sets the amount of memory to use for the Spark
driver process, which is the process running the user code that creates
RDDs, accumulators, and broadcast variables, and manages their storage
in memory. Here it’s being set to 4 gigabytes.

**Executor memory**: This sets the memory size for Spark executors.
Executors are worker nodes’ processes in charge of running the tasks in
a Spark job. Here each executor is allocated 50 gigabytes of memory.

**offHeap size**: This enables off-heap memory allocation. This setting
is false by default and uses on-heap memory. Off-heap memory can be used
for certain operations like shuffle and can help mitigate garbage
collection overhead.

.. ipython:: python

    from pyspark.sql import SparkSession
    from pyspark import SparkContext
    
    from pyspark.sql import SparkSession
    
    spark = SparkSession.builder \
        .appName("StructureAnalysis") \
        .config("spark.driver.memory", "4g") \
        .config("spark.executor.memory", "50g") \
        .config("spark.memory.offHeap.enabled", True) \
        .config("spark.memory.offHeap.size","8g") \
        .getOrCreate()
    
    sc = spark.sparkContext


We can then read in the sequence file through the spark context. With
large computing resources, the full sequence file can be analyzed quite
efficiently. ‘input/full/’ is the path to the folder that contains the
partitioned sequence files. We will sample 10% of the entire PDB
archive.

.. ipython:: python

    mmtf_full = sc.sequenceFile('input/full/')
    mmtf_sample = mmtf_full.sample(False, 0.1).persist()

**mmtf_full** and **mmtf_sample** are RDDs (Resilient Distributed
Dataset) on which you can apply PySpark functions such as map, filter,
or count. Each entry is a a tuple of the form **(‘PDB_ID’,
byte_array)**. For example:

.. ipython:: python

    sample = mmtf_sample.take(1)
    sample
                                                                                  

.. parsed-literal::

    [('1GQ8',
      bytearray(b'\x1f\x8b...'))]


:func:`.parseMMTF` can parse MMTF files, structure objects, and even byte
arrays. In order to parse the byte-array to a ProDy Atom Group:

.. ipython:: python

    struc_byte_array = sample[0][1]
    
    structure = parseMMTF(struc_byte_array)
    showProtein(structure)


.. parsed-literal::

    @> 2770 atoms and 1 coordinate set(s) were parsed in 0.01s.

.. parsed-literal::

    <py3Dmol.view at 0x7fb4d4671190>



We can also write a function that can be mapped to the entire dataset to
retrieve information, whether to extract structures, or to perform
exploratory data analysis. Here, we will go through the entire sampled
dataset and return the top 10 PDB structures sorted by number of atoms.

.. ipython:: python

    def get_num_atoms(tuple_entry):
        structure = parseMMTF(tuple_entry[1])
        
        return (structure.numAtoms(),structure.getTitle())
    
    top_10_entries = mmtf_sample.map(get_num_atoms) 
    top_10_entries.takeOrdered(10, key=lambda x: -x[0])


.. parsed-literal::

    [(1324599, '8G3D'),
     (1265783, '7UNG'),
     (1234795, '5Y6P'),
     (1218816, '7SQC'),
     (1013398, '6KGX'),
     (662907, '7AS5'),
     (617133, '7R5J'),
     (590573, '4V9P'),
     (511395, '4V8P'),
     (411205, '4U3U')]



We can now terminate the spark session.

.. ipython:: python

    sc.stop()
