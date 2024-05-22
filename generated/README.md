# Research Software Directory

This repo contains [data.yml](data.yml) which is the data used by
- [SSC Research Software Directory](https://www.ssc.uni-heidelberg.de/en/research-software-directory)
- [SSC Research Software Directory Visualisation](https://ssciwr.github.io/research-software-directory-visualization)

Whenever [data.yml](data.yml) is updated, the HTML code below is automatically updated, which can be copy&pasted (click icon in top right of code block below) into drupal:

*Note: this README.md file is automatically generated from the [README.md.j2](README.md.j2) jinja template - edit the template if you want to modify this README!*

```html
<p>There is a lot of research software that is being developed and used at Heidelberg University. This Research Software Directory is meant to be a comprehensive but by no means complete list of such software, and may aid you in identifying collaboration partners.</p>
<table>
  <thead>
    <tr>
      <th>Software</th>
      <th>Field of Study</th>
      <th>Research Group</th>
      <th>DOI</th>
      <th>Link</th>
    </tr>
  </thead>
  <tbody>
    <tr>
      <td>
        <p><strong>A-SLOTH</strong></p>
        <p>The semi-analytical model A-SLOTH (Ancient Stars and Local Observables by Tracing Halos) is the first public code that connects the formation of the first stars and galaxies to observables. The model is based on dark matter merger trees that can either be generated based on Extended Press-Schechter theory or that can be imported from dark matter simulations.&nbsp;</p>
      </td>
      <td>Physics, Astrophysics</td>
      <td>Klessen group</td>
      <td>
        <p>10.21105/joss.04417</p>
        <p>10.3847/1538-4357/ac7150</p>
      </td>
      <td><a href="https://gitlab.com/thartwig/asloth">https://gitlab.com/thartwig/asloth</a></td>
    </tr>
    <tr>
      <td>
        <p><strong>ACRONYM</strong></p>
        <p>ACRONYM is an explicit, highly-flexible particle-in-cell code for simulations of 1d, 2d, or 3d spatial domains with a number of charged particle species. Particle motion is relativistically correct and computed using the Boris push. Temporal evolution of the fields can be computed for electrostatic or electromagnetic plasma models using a number of different field solvers. Interpolation between field- and particle-quantities can be done through a choice of plugable shape functions.</p>
      </td>
      <td>Physics, Astrophysics</td>
      <td>Spanier group</td>
      <td>
        <p>10.1007/978-3-642-23869-7_1</p>
      </td>
      <td><a href="https://www.ita.uni-heidelberg.de/~fspanier/software/acronym.shtml?lang…">https://www.ita.uni-heidelberg.de/~fspanier/software/acronym.shtml?lang…</a></td>
    </tr>
    <tr>
      <td>
        <p><strong>AFTERLIVE</strong></p>
        <p>Afterlive (A Fourier-based Tool in the Electrostatic limit for the Rapid Low-noise Integration of the Vlasov Equation) allows for a fast simulation of electromagnetic plasmas. It uses a hybrid method between a gridded Eulerian description and Lagrangian meta-particles, reconstructing the distribution function f in every time step for each species. This interpolation method combines meta-particles with different weights in such a way that particles with large weight do not drown out particles that represent small contributions to the phase space density. These core properties allow the use of a much larger range of macro factors and can thus represent a much larger dynamic range in phase space density.</p>
      </td>
      <td>Physics, Astrophysics</td>
      <td>Spanier group</td>
      <td>
        <p>10.1016/j.cpc.2018.04.014</p>
        <p>10.17632/nt79b5nfsc.1</p>
      </td>
      <td><a href="https://www.ita.uni-heidelberg.de/~fspanier/software/afterlive.shtml?la…">https://www.ita.uni-heidelberg.de/~fspanier/software/afterlive.shtml?la…</a></td>
    </tr>
    <tr>
      <td>
        <p><strong>AMMICO</strong></p>
        <p>AMMICO, the AI Media and Misinformation Content Analysis Tool, extracts data from images such as social media posts that contain an image and text. The analysis can generate a very large number of features, for example: Text extraction from the images, language detection, translation into English or other languages, sentiment analysis, named entity recognition, image caption generation based on the image content, image matching, visual question answering, face recognition, face mask detection, emotion recognition, and more.</p>
      </td>
      <td>Social Sciences</td>
      <td>Dumistrescu group/SSC</td>
      <td>
        <p>10.31235/osf.io/v8txj</p>
      </td>
      <td><a href="https://github.com/ssciwr/AMMICO">https://github.com/ssciwr/AMMICO</a></td>
    </tr>
    <tr>
      <td>
        <p><strong>ArepoVTK</strong></p>
        <p>ArepoVTK is a visualization library designed to produce high quality, presentation-ready visualizations of hydrodynamic simulations run with the AREPO code. It can render single images as well as frame sequences for movies from cosmological simulations such as IllustrisTNG as well as idealized test problems and other computational fluid dynamics simulations.</p>
      </td>
      <td>Physics, Astrophysics, Visual Computing</td>
      <td>Nelson group</td>
      <td>
        <p>10.1093/mnras/sts595</p>
      </td>
      <td><a href="https://github.com/nelson-group/ArepoVTK">https://github.com/nelson-group/ArepoVTK</a></td>
    </tr>
    <tr>
      <td>
        <p><strong>CARNIVAL</strong></p>
        <p>CARNIVAL (CAusal Reasoning for Network identification using Integer VALue programming) is a method for the identification of upstream regulatory signalling pathways from downstream gene expression (GEX).The aim of the CARNIVAL is to identify a subset of interactions from a prior knowledge network that represent potential regulated pathways linking known or potential targets of perturbation towards active transcription factors derived from GEX data.&nbsp;</p>
      </td>
      <td>Biology, Bioinformatics</td>
      <td>Saez-Rodriguez group</td>
      <td>
        <p>10.1038/s41540-019-0118-z</p>
      </td>
      <td><a href="https://saezlab.github.io/CARNIVAL/">https://saezlab.github.io/CARNIVAL/</a></td>
    </tr>
    <tr>
      <td>
        <p><strong>CONAN</strong></p>
        <p>Contact analysis can be extremely useful to understand structure and dynamics of molecules and molecular complexes as they are atomistic-resolved descriptions of the interactions occurring within and between the investigated molecules. CONAN is a tool developed for the statistical and dynamical analysis of contacts along molecular dynamics trajectories. It uses a combination of open-source tools for the computation of contacts along trajectories and for the creation of images and videos.&nbsp;</p>
      </td>
      <td>Physical Chemistry, Biochemistry, Biophysics</td>
      <td>Gräter group</td>
      <td>
        <p></p>
      </td>
      <td><a href="https://hits-mbm.github.io/conan/">https://hits-mbm.github.io/conan/</a></td>
    </tr>
    <tr>
      <td>
        <p><strong>ConDec</strong></p>
        <p>ConDec supports developers in continuously managing decision knowledge. To reduce intrusiveness and additional effort, ConDec integrates into the workflows and tools of agile software development. ConDec enables lightweight capturing and use of decision knowledge from various artifacts and reduces the developers' effort through automatic text classification, recommendation, and nudging mechanisms. To enable access and use of distributed decision knowledge documentation, ConDec builds up a knowledge graph of requirements, decision knowledge, code, and other artifacts, offering interactive knowledge views. To support high quality, ConDec offers the rationale backlog, the definition of done for knowledge documentation, the knowledge dashboard and metrics for intra-rationale completeness and decision coverage of requirements and code. ConDec supports consistent changes through change impact analysis.</p>
      </td>
      <td>Computer Science, Software Engineering</td>
      <td>Paech group</td>
      <td>
        <p>10.5281/zenodo.7953646</p>
        <p>10.5281/zenodo.7948296</p>
        <p>10.5281/zenodo.6866299</p>
        <p>10.5281/zenodo.7958350</p>
        <p>10.5281/zenodo.7954331</p>
        <p>10.5281/zenodo.7955051</p>
      </td>
      <td><a href="https://github.com/cures-hub">https://github.com/cures-hub</a></td>
    </tr>
    <tr>
      <td>
        <p><strong>COPASI</strong></p>
        <p>COPASI is a software application for simulation and analysis of biochemical networks and their dynamics. COPASI supports models in the SBML standard and can simulate their behavior using ODEs, SDEs, or Gillespie's stochastic simulation algorithm. Arbitrary discrete events can be included in such simulations. COPASI provides a set analysis methods and parameter estimation.</p>
      </td>
      <td>Biology, Bioinformatics</td>
      <td>Kummer group</td>
      <td>
        <p>10.1093/bioinformatics/btl485</p>
      </td>
      <td><a href="http://copasi.org/">http://copasi.org/</a></td>
    </tr>
    <tr>
      <td>
        <p><strong>cosmoRate</strong></p>
        <p>cosmoRate is a&nbsp;semi-analytic code for the merger rate density evolution of binary compact object across cosmic time.</p>
      </td>
      <td>Physics, Astrophysics</td>
      <td>Mapelli group</td>
      <td>
        <p>10.1093/mnras/stad1860</p>
        <p>10.1093/mnras/stab280</p>
      </td>
      <td><a href="https://gitlab.com/Filippo.santoliquido/cosmo_rate_public">https://gitlab.com/Filippo.santoliquido/cosmo_rate_public</a></td>
    </tr>
    <tr>
      <td>
        <p><strong>dantro</strong></p>
        <p>dantro is a Python package that provides a uniform interface for hierarchically structured and semantically heterogeneous data. It provides an integrated, automated and configurable data processing pipeline from data handling over arbitrary transformation to publication quality data visualization.</p>
      </td>
      <td>Data Science</td>
      <td>Roth group</td>
      <td>
        <p>10.21105/joss.02316</p>
      </td>
      <td><a href="https://github.com/utopia-foss/dantro">https://github.com/utopia-foss/dantro</a></td>
    </tr>
    <tr>
      <td>
        <p><strong>decoupleR</strong></p>
        <p>There are many methods that allow us to extract biological activities from omics data.&nbsp;decoupleR&nbsp;is a Bioconductor package containing different statistical methods to extract biological signatures from prior knowledge within a unified framework. Additionally, it incorporates methods that take into account the sign and weight of network interactions.&nbsp;decoupleR&nbsp;can be used with any omic, as long as its features can be linked to a biological process based on prior knowledge. For example, in transcriptomics gene sets regulated by a transcription factor, or in phospho-proteomics phosphosites that are targeted by a kinase. There are both R and python versions.&nbsp;</p>
      </td>
      <td>Biology, Bioinformatics</td>
      <td>Saez-Rodriguez group</td>
      <td>
        <p>10.1093/bioadv/vbac016</p>
      </td>
      <td><a href="https://saezlab.github.io/decoupleR/">https://saezlab.github.io/decoupleR/</a></td>
    </tr>
    <tr>
      <td>
        <p><strong>deal.II</strong></p>
        <p>deal.II is an open source finite element library supporting the creation of finite element codes and an open community of users and developers. deal.II enables rapid development through taking care of the details of grid handling and refinement, handling of degrees of freedom, input of meshes and output of results in graphics formats, among others.</p>
      </td>
      <td>Computer Science, Applied Mathematics</td>
      <td>Kanschat group</td>
      <td>
        <p>10.1515/jnma-2023-0089</p>
        <p>10.1016/j.camwa.2020.02.022</p>
      </td>
      <td><a href="https://www.dealii.org/">https://www.dealii.org/</a></td>
    </tr>
    <tr>
      <td>
        <p><strong>DUNE</strong></p>
        <p>DUNE, the Distributed and Unified Numerics Environment is a modular toolbox for solving partial differential equations&nbsp;with grid-based methods. It supports the easy implementation of methods like Finite Elements, Finite Volumes, and also Finite Differences.</p>
      </td>
      <td>Computer Science, Applied Mathematics</td>
      <td>Bastian group</td>
      <td>
        <p>10.1016/j.camwa.2020.06.007</p>
      </td>
      <td><a href="https://www.dune-project.org/">https://www.dune-project.org/</a></td>
    </tr>
    <tr>
      <td>
        <p><strong>FASTCLUSTER</strong></p>
        <p>FASTCLUSTER is an open-source population-synthesis code for binary black hole dynamics in star clusters. You can use FASTCLUSTER to study the hierarchical formation of binary black holes.</p>
      </td>
      <td>Physics, Astrophysics</td>
      <td>Mapelli group</td>
      <td>
        <p>10.1093/mnras/stac422</p>
        <p>10.1093/mnras/stab1334</p>
      </td>
      <td><a href="https://gitlab.com/micmap/fastcluster_open">https://gitlab.com/micmap/fastcluster_open</a></td>
    </tr>
    <tr>
      <td>
        <p><strong>FoNo</strong></p>
        <p>FolioNumbering (FoNo) is a Python tool to rename image files in a meaningful way.</p>
      </td>
      <td>Text Edition, Corpus Linguistics, Historical Lexicology, Linked Open Data</td>
      <td>Tittel group/ALMA project</td>
      <td>
        <p></p>
      </td>
      <td><a href="https://gitlab.hadw-bw.de/alma/alma/-/tree/main/fono">https://gitlab.hadw-bw.de/alma/alma/-/tree/main/fono</a></td>
    </tr>
    <tr>
      <td>
        <p><strong>Force Distribution Analysis</strong></p>
        <p>Force distribution analysis (FDA) is a method to detect and follow force and stress propagation in proteins, reminiscent of Finite Element Analysis used to engineer macroscopic structures. The method is based on Molecular Dynamics simulations during which we directly calculate forces between each atom pair in a system.&nbsp;</p>
      </td>
      <td>Physical Chemistry, Biochemistry, Biophysics</td>
      <td>Gräter group</td>
      <td>
        <p>10.1186/2046-1682-6-5</p>
      </td>
      <td><a href="https://github.com/HITS-MBM/gromacs-fda">https://github.com/HITS-MBM/gromacs-fda</a></td>
    </tr>
    <tr>
      <td>
        <p><strong>FISHFactor</strong></p>
        <p>FISHFactor is a&nbsp;probabilistic factor model for spatial transcriptomics data with subcellular resolution.&nbsp;It provides a&nbsp;non-negative, spatially informed factor analysis model with a Poisson point process likelihood to model single-molecule resolved data, as obtained for example from multiplexed fluorescence in-situ hybridization methods, and allows to integrate multiple cells by jointly inferring cell-specific factors and a weight matrix that is shared between cells.</p>
      </td>
      <td>Biology, Bioinformatics</td>
      <td>Velten group</td>
      <td>
        <p>10.1093/bioinformatics/btad183</p>
      </td>
      <td><a href="https://github.com/bioFAM/FISHFactor">https://github.com/bioFAM/FISHFactor</a></td>
    </tr>
    <tr>
      <td>
        <p><strong>GASOLINE</strong></p>
        <p>Particle Hydrodynamics Have Never Been Smoother: Gasoline is a modern TREESPH code for solving the equations of gravity and hydrodynamics in astrophysical problems.</p>
      </td>
      <td>Physics, Astrophysics</td>
      <td>Buck group</td>
      <td>
        <p></p>
      </td>
      <td><a href="https://gasoline-code.com/">https://gasoline-code.com/</a></td>
    </tr>
    <tr>
      <td>
        <p><strong>GaussianVariationalInference.jl</strong></p>
        <p>GaussianVariationalInference.jl is a Julia package for approximating a posterior distribution with a full-covariance Gaussian distribution by optimising a variational lower bound. The main focus of this package is to provide a method for approximating a target posterior with a Gaussian that does not need tuning learning rates (step sizes) and converges reliably.</p>
      </td>
      <td>Machine Learning</td>
      <td>Astroinformatics (HITS)</td>
      <td>
        <p>10.1007/s10044-015-0496-9</p>
      </td>
      <td><a href="https://github.com/ngiann/GaussianVariationalInference.jl">https://github.com/ngiann/GaussianVariationalInference.jl</a></td>
    </tr>
    <tr>
      <td>
        <p><strong>GPCC</strong></p>
        <p>A Gaussian process cross-correlation approach to time delay estimation for reverberation mapping of active galactic nuclei. GPCC is a probabilistic alternative to the Interpolated Cross Correlation function (ICCF). Advantages over the ICCF include generation of a probability distribution for the delay, possibility to incorporate a prior on the delay, predictions for out-of-sample data.</p>
      </td>
      <td>Astrophysics</td>
      <td>Astroinformatics (HITS)</td>
      <td>
        <p>10.1051/0004-6361/202345932</p>
        <p>ascl code record ascl:2303.006</p>
      </td>
      <td><a href="https://github.com/HITS-AIN/GPCC.jl">https://github.com/HITS-AIN/GPCC.jl</a></td>
    </tr>
    <tr>
      <td>
        <p><strong>GROmaps</strong></p>
        <p>GROmaρs is an open-source GROMACS-based toolset for the calculation and comparison of density maps from molecular dynamics simulations.</p>
      </td>
      <td>Physical Chemistry, Biochemistry, Biophysics</td>
      <td>Gräter group</td>
      <td>
        <p>10.1016/j.bpj.2018.11.3126</p>
      </td>
      <td><a href="https://mptg-cbp.github.io/gromaps.html">https://mptg-cbp.github.io/gromaps.html</a></td>
    </tr>
    <tr>
      <td>
        <p><strong>gsf</strong></p>
        <p>gsf applies Gaussian Mixture Models in the stellar kinematic space of normalized angular momentum and binding energy on NIHAO high resolution galaxies to separate the stars into multiple components. The gsf analysis package assumes that the simulation snapshot has been pre-processed with a halo finder. It is based on pynbody and the scikit-learn Python package for Machine Learning; after loading, orienting, and transforming a simulation snapshot to physical units, it runs the clustering algorithm and performs the direct N-body gravity force using all the particles in the given halo.</p>
      </td>
      <td>Physics, Astrophysics</td>
      <td>Buck group</td>
      <td>
        <p>10.1093/mnras/sty1022</p>
      </td>
      <td><a href="https://github.com/aobr/gsf">https://github.com/aobr/gsf</a></td>
    </tr>
    <tr>
      <td>
        <p><strong>HELIOS++</strong></p>
        <p>HELIOS++ is a general-purpose software package for simulation of terrestrial, mobile and airborne laser scanning surveys. Virtual laser scanning is a tool to create simulated point cloud data, as would be acquired by a LiDAR sensor. Such data may be used to complement real data, where data acquisition is not feasible due to economical or logistic constraints or where it is impossible, e.g. when simulating a sensor that does not exist. HELIOS++ allows the simulation of laser scanning on different platforms (airborne, UAV-based, terrestrial mobile and static) and using different data types to represent the 3D scene, including triangular meshes, digital elevation rasters, voxel grids and point clouds.</p>
      </td>
      <td>Geoscience, Geoinformatics, Geodesy</td>
      <td>Höfle group</td>
      <td>
        <p>10.1016/j.rse.2021.112772</p>
        <p>10.1109/ACCESS.2022.3211072</p>
      </td>
      <td><a href="https://github.com/3dgeo-heidelberg/helios">https://github.com/3dgeo-heidelberg/helios</a></td>
    </tr>
    <tr>
      <td>
        <p><strong>Hiflow3</strong></p>
        <p>HiFlow3 is a multi-purpose finite element software providing powerful tools for efficient and accurate solution of a wide range of problems modeled by partial differential equations. Based on object-oriented concepts and the full capabilities of C++ the HiFlow³ project follows a modular and generic approach for building efficient parallel numerical solvers.</p>
      </td>
      <td>Computer Science, Applied Mathematics</td>
      <td>Heuveline group</td>
      <td>
        <p>10.11588/emclpp.2017.06.42879</p>
      </td>
      <td><a href="http://hiflow3.org/">http://hiflow3.org/</a></td>
    </tr>
    <tr>
      <td>
        <p><strong>HyperHDG</strong></p>
        <p>HyperHDG is a C++ based library implementing hybrid discontinuous Galerkin methods on extremely general domains: that is, all standard (volume) domains, graphs, networks of surfaces, and several other types of "domains" that can be interpreted as hypergraphs.</p>
      </td>
      <td>Computer Science, Applied Mathematics</td>
      <td>Kanschat group</td>
      <td>
        <p>10.1051/m2an/2022011</p>
      </td>
      <td><a href="https://github.com/HyperHDG/HyperHDG">https://github.com/HyperHDG/HyperHDG</a></td>
    </tr>
    <tr>
      <td>
        <p><strong>ilastik</strong></p>
        <p>ilastik the interactive learning and segmentation toolkit. Leverage machine learning algorithms to easily segment, classify, track and count your cells or other experimental data. Most operations are interactive, even on large datasets: you just draw the labels and immediately see the result.</p>
      </td>
      <td>Computer Vision, Machine Learning</td>
      <td>Kreshuk group/Hamprecht group</td>
      <td>
        <p>10.1109/ISBI.2011.5872394</p>
        <p>10.1007/978-3-319-28549-8_8</p>
      </td>
      <td><a href="https://www.ilastik.org/">https://www.ilastik.org/</a></td>
    </tr>
    <tr>
      <td>
        <p><strong>KBbox</strong></p>
        <p>KBbox (A Toolbox of Computational Methods for Studying the Kinetics of Molecular Binding) helps researchers decide which method for studying the kinetics of binding processes and predicting their associated rate constants might be suitable for their projects. KBbox is a web server that guides users in choosing the methods they should consider on the basis of the information they wish to obtain, the data they currently have available, and the computational resources to which they have access. KBbox provides information on the toolbox of available methods, their associated software tools, an expanding list of curated examples of published applications, and tutorials explaining how to apply some of the methods. It has been designed to allow the easy addition of new methods, tools, and examples as they are developed and published.&nbsp;</p>
      </td>
      <td>Biology, Bioinformatics</td>
      <td>Wade group</td>
      <td>
        <p>10.1021/acs.jcim.9b00485</p>
      </td>
      <td><a href="https://kbbox.h-its.org/toolbox/">https://kbbox.h-its.org/toolbox/</a></td>
    </tr>
    <tr>
      <td>
        <p><strong>KIMMDY</strong></p>
        <p>KIMMDY (Kinetic Monte Carlo / Molecular Dynamics) enables covalent bond breakages in all-atom Molecular Dynamics (MD) simulations. The bond rupture rates are calculated based on the interatomic distances in the MD simulation and then serve as an input for a Kinetic Monte Carlo step. This hybrid approach allows to study breakages in large molecular systems bridging various time scales between MD and the rupture processes.&nbsp;</p>
      </td>
      <td>Physical Chemistry, Biochemistry, Biophysics</td>
      <td>Gräter group</td>
      <td>
        <p>10.1021/acs.jctc.9b00786</p>
      </td>
      <td><a href="https://github.com/HITS-MBM/kimmdy">https://github.com/HITS-MBM/kimmdy</a></td>
    </tr>
    <tr>
      <td>
        <p><strong>L-RIP / RIPLig</strong></p>
        <p>Rotamerically Induced Perturbations with Langevin dynamics (L-RIP) and Rotamerically Induced Perturbations of a pseudo-Ligand (RIPlig) are two non-equilibrium MD approaches aimed at an identification of slow conformational changes of a binding site. They are an extention of the Rotamerically Induced Perturbation (RIP) MD approach, which employs perturbation of side-chain torsional motion for initiating large-scale protein movement. RIP enables protein flexibility that normally occurs on the microsecond or longer time scale (for example, distortion of secondary and tertiary structure) to be explored. In Langevin-RIP (L-RIP) approach perturbation is applied to the side chains of the binding site residues, but the original RIP procedure is modified by using a Langevin thermostat in the simulations instead of the original constant-energy conditions. This ensures that the average protein temperature is preserved during the MD simulations. Additionally, the friction term enables the elements of the protein structure with fast relaxation times to equilibrate within each perturbation pulse. This makes it possible to increase the number of pulses and, thus, enhance sampling of slow conformational variations. In RIPlig, the perturbation procedure is applied not to each protein residue, but to an artificial molecule placed in the binding pocket with repetitions starting with different positions of the molecule, which enables flexible regions of the binding site to be identified in a small number of MD simulations of several picoseconds duration.</p>
      </td>
      <td>Biology, Bioinformatics</td>
      <td>Wade group</td>
      <td>
        <p>10.1021/acs.jctc.6b00101</p>
      </td>
      <td><a href="https://mcm.h-its.org/lrip-riplig/doc/main.html">https://mcm.h-its.org/lrip-riplig/doc/main.html</a></td>
    </tr>
    <tr>
      <td>
        <p><strong>LIANA+</strong></p>
        <p>LIANA+ is an efficient framework that integrates and extends existing methods and knowledge to study cell-cell communication in single-cell, spatially-resolved, and multi-modal omics data.</p>
      </td>
      <td>Biology, Bioinformatics</td>
      <td>Saez-Rodriguez group</td>
      <td>
        <p>10.1038/s41467-022-30755-0</p>
      </td>
      <td><a href="https://liana-py.readthedocs.io/en/latest/">https://liana-py.readthedocs.io/en/latest/</a></td>
    </tr>
    <tr>
      <td>
        <p><strong>Luminescence</strong></p>
        <p>The Luminescence package is a collection of various R functions for the purpose of Luminescence dating data analysis. This includes, amongst others, data import, export, application of age models, curve deconvolution, sequence analysis and plotting of equivalent dose distributions.</p>
      </td>
      <td>Geoscience, Geoinformatics, Applied Physics</td>
      <td>Kreutzer group</td>
      <td>
        <p>Ancient TL, 30(1), 1-8 (2012)</p>
        <p>10.5281/zenodo.596252</p>
      </td>
      <td><a href="https://cran.r-project.org/web/packages/Luminescence/index.html">https://cran.r-project.org/web/packages/Luminescence/index.html</a></td>
    </tr>
    <tr>
      <td>
        <p><strong>MapSwipe</strong></p>
        <p>MapSwipe harnesses the collective strength of volunteers to actively contribute to geospatial data projects. From identifying infrastructure to tracking environmental changes and validating map data, MapSwipers help improve map data across the world.</p>
      </td>
      <td>Geoscience, Geoinformatics</td>
      <td>Zipf group</td>
      <td>
        <p></p>
      </td>
      <td><a href="https://mapswipe.org/en/">https://mapswipe.org/en/</a></td>
    </tr>
    <tr>
      <td>
        <p><strong>MCTDH</strong></p>
        <p>Multi Configuration Time Dependent Hartree (MCTDH) is an algorithm to solve the time-dependent Schrödinger equation for multidimensional dynamical systems consisting of distinguishable particles, such as the quantal motion of the nuclei of a molecular system evolving on one or several coupled electronic potential energy surfaces.</p>
      </td>
      <td>Theoretical chemistry, Physical Chemistry</td>
      <td>Vendrell group</td>
      <td>
        <p>10.1016/0009-2614(90)87014-I</p>
        <p>10.1016/S0370-1573(99)00047-2</p>
      </td>
      <td><a href="https://www.pci.uni-heidelberg.de/cms/mctdh.html">https://www.pci.uni-heidelberg.de/cms/mctdh.html</a></td>
    </tr>
    <tr>
      <td>
        <p><strong>MD-IFP</strong></p>
        <p>MD-IFP provides an MD trajectory analysis using protein-ligand Interaction Fingerprints. It is a Python Workflow for the Generation and Analysis of Protein-Ligand Interaction Fingerprints from Molecular Dynamics trajectories.</p>
      </td>
      <td>Biology, Bioinformatics</td>
      <td>Wade group</td>
      <td>
        <p>10.1063/5.0019088</p>
      </td>
      <td><a href="https://github.com/HITS-MCM/MD-IFP">https://github.com/HITS-MCM/MD-IFP</a></td>
    </tr>
    <tr>
      <td>
        <p><strong>MEFISTO</strong></p>
        <p>MEFISTO provides an unsupervised approach to integrate multi-modal data with continuous structures among the samples, e.g. given by spatial or temporal relationships. The aim of MEFISTO is to exploit such relationships between samples in the dimensionality reduction and disentangle smooth sources of variation given by factors that change gradually along the covariate and other source of variation that are independent of the covariate. Furthermore, it enables to interpolate/extrapolate to unseen timepoints or locations.</p>
      </td>
      <td>Biology, Bioinformatics</td>
      <td>Velten group</td>
      <td>
        <p>10.1038/s41592-021-01343-9</p>
      </td>
      <td><a href="https://biofam.github.io/MOFA2/MEFISTO.html">https://biofam.github.io/MOFA2/MEFISTO.html</a></td>
    </tr>
    <tr>
      <td>
        <p><strong>MISTy</strong></p>
        <p>The advancement of technologies for measurement of highly multiplexed spatial data require the development of scalable methods that can leverage the availability of the spatial context. Multiview Intercellular SpaTial modeling framework (MISTy) is an explainable machine learning framework for knowledge extraction and analysis of single-cell, highly multiplexed, spatially resolved data.</p>
      </td>
      <td>Biology, Bioinformatics</td>
      <td>Saez-Rodriguez group</td>
      <td>
        <p>10.1186/s13059-022-02663-5</p>
      </td>
      <td><a href="https://saezlab.github.io/mistyR/">https://saezlab.github.io/mistyR/</a></td>
    </tr>
    <tr>
      <td>
        <p><strong>MOBSE</strong></p>
        <p>MOBSE investigates the demography of merging BHBs. A customized version of the binary stellar evolution code BSE (ascl:1303.014), MOBSE includes metallicity-dependent prescriptions for mass-loss of massive hot stars and upgrades for the evolution of single and binary massive stars.</p>
      </td>
      <td>Physics, Astrophysics</td>
      <td>Mapelli group</td>
      <td>
        <p>10.1093/mnras/stx2123</p>
        <p>ads bibcode 2023ascl.soft06010G</p>
      </td>
      <td><a href="https://gitlab.com/micmap/mobse_open">https://gitlab.com/micmap/mobse_open</a></td>
    </tr>
    <tr>
      <td>
        <p><strong>Molsurfer</strong></p>
        <p>MolSurfer has been designed to assist the analysis of the structures and physico-chemical properties of macromolecular interfaces. MolSurfer provides a coupled display of two-dimensional (2D) maps of the interfaces generated with the ADS software and a three-dimensional (3D) view of the macromolecular structure in the Java PDB viewer, WebMol. The interfaces are analytically defined and properties such as electrostatic potential or hydrophobicity are projected on to them.&nbsp;</p>
      </td>
      <td>Biology, Bioinformatics</td>
      <td>Wade group</td>
      <td>
        <p>10.1093/nar/gkg588</p>
        <p>10.1016/s0968-0004(99)01412-7</p>
      </td>
      <td><a href="https://molsurfer.h-its.org/">https://molsurfer.h-its.org/</a></td>
    </tr>
    <tr>
      <td>
        <p><strong>MOFA</strong></p>
        <p>Multi-omics factor analysis (MOFA)&nbsp;provides a general framework for the integration of multi-omic data sets in an unsupervised fashion. Intuitively, MOFA can be viewed as a versatile and statistically rigorous generalization of principal component analysis to multi-omics data.&nbsp;</p>
      </td>
      <td>Biology, Bioinformatics</td>
      <td>Velten group</td>
      <td>
        <p>10.15252/msb.20178124</p>
        <p>10.1186/s13059-020-02015-1</p>
      </td>
      <td><a href="https://biofam.github.io/MOFA2/index.html">https://biofam.github.io/MOFA2/index.html</a></td>
    </tr>
    <tr>
      <td>
        <p><strong>ohsome</strong></p>
        <p>ohsome is an OpenStreetMap History Data Analytics Platform. Our aim is to make OpenStreetMap’s full-history data more easily accessible for various kinds of data analytics tasks on a global scale. Applications of the ohsome platform range from web dashboards over data quality assessment to custom data analysis.</p>
      </td>
      <td>Geoscience, Geoinformatics</td>
      <td>Zipf group</td>
      <td>
        <p></p>
      </td>
      <td><a href="https://heigit.org/big-spatial-data-analytics-en/ohsome/">https://heigit.org/big-spatial-data-analytics-en/ohsome/</a></td>
    </tr>
    <tr>
      <td>
        <p><strong>Openrouteservice</strong></p>
        <p>The openrouteservice API provides global spatial services by consuming user-generated and collaboratively collected free geographic data directly from OpenStreetMap. The following services are available: Directions - returns a route between two or more locations for a selected profile with customizable additional settings and instructions; Isochrones - obtains areas of reachability from given locations; Matrix - computes one-to-many, many-to-one or many-to-many routes for any mode of transport provided by openrouteservice.</p>
      </td>
      <td>Geoscience, Geoinformatics</td>
      <td>Zipf group</td>
      <td>
        <p></p>
      </td>
      <td><a href="https://openrouteservice.org/">https://openrouteservice.org/</a></td>
    </tr>
    <tr>
      <td>
        <p><strong>PICPANTHER</strong></p>
        <p>PICPANTHER implements the relativistic moment implicit particle-in-cell method. The implicit electric field equation is solved using the GMRES algorithm with operators represented as sparse matrices. For each particle, the implicit equation of motion is solved via a robust Newton-Krylov scheme. Parallelization is achieved using simple domain decomposition, resulting in good scalability. PICPANTHER makes use of advanced numerical techniques (GMRES, Newton-Krylov) to implicitly solve relativistic versions of the movement and field equations of a PiC code.</p>
      </td>
      <td>Physics, Astrophysics</td>
      <td>Spanier group</td>
      <td>
        <p>10.1016/j.cpc.2014.11.010,</p>
        <p>10.17632/gcy5p3ddfb.1</p>
      </td>
      <td><a href="https://www.ita.uni-heidelberg.de/~fspanier/software/picpanther.shtml?l…">https://www.ita.uni-heidelberg.de/~fspanier/software/picpanther.shtml?l…</a></td>
    </tr>
    <tr>
      <td>
        <p><strong>PIPSA / webPIPSA</strong></p>
        <p>PIPSA allows a comparison of the electrostatic interaction properties of proteins. It permits the classification of proteins according to their interaction properties. PIPSA may assist in function assignment, the estimation of binding properties and enzyme kinetic parameters.&nbsp;</p>
      </td>
      <td>Biology, Bioinformatics</td>
      <td>Wade group</td>
      <td>
        <p>10.1002/(sici)1097-0134(19991115)37:3&lt;379::aid-prot6&gt;3.0.co;2-k</p>
        <p>10.1002/qua.1204</p>
        <p>10.1186/1471-2105-8-373</p>
        <p>10.1093/nar/gkn181</p>
      </td>
      <td><a href="https://pipsa.h-its.org/pipsa/">https://pipsa.h-its.org/pipsa/</a></td>
    </tr>
    <tr>
      <td>
        <p><strong>PlantSeg</strong></p>
        <p>PlantSeg is a tool for cell instance aware segmentation in densely packed 3D volumetric images. The pipeline uses a two stages segmentation strategy (Neural Network + Segmentation). The pipeline is tuned for plant cell tissue acquired with confocal and light sheet microscopy. Pre-trained models are provided.</p>
      </td>
      <td>Image Analysis, Machine Learning, Biology</td>
      <td>Hamprecht group</td>
      <td>
        <p>10.7554/eLife.57613</p>
      </td>
      <td><a href="https://github.com/hci-unihd/plant-seg">https://github.com/hci-unihd/plant-seg</a></td>
    </tr>
    <tr>
      <td>
        <p><strong>POLARIS</strong></p>
        <p>POLARIS (POLArized RadIation Simulator) is a 3D Monte-Carlo radiative transfer code allows to simulate intensity and polarization of light emerging from analytical astrophysical models as well as complex magneto-hydrodynamic simulations on various grids. POLARIS is capable to perform dust heating, -emission, -scattering, -grain alignment, line radiative transfer, and synchrotron simulations.</p>
      </td>
      <td>Physics, Astrophysics</td>
      <td>Klessen group</td>
      <td>
        <p>10.1051/0004-6361/201424930</p>
        <p>10.1051/0004-6361/201629001</p>
      </td>
      <td><a href="https://portia.astrophysik.uni-kiel.de/polaris/">https://portia.astrophysik.uni-kiel.de/polaris/</a></td>
    </tr>
    <tr>
      <td>
        <p><strong>py4dgeo</strong></p>
        <p>py4dgeo is a C++ library with Python bindings for change analysis in multitemporal and 4D point clouds. Topographic 3D/4D point clouds are omnipresent in geosciences, environmental, ecological and archaeological sciences, robotics, and many more fields and applications. Technology to capture such data using laser scanning and photogrammetric techniques have evolved into standard tools. Dense time series of topographic point clouds are becoming increasing available and require tools for automatic analysis. Moreover, methods considering the full 4D (3D space + time) data are being developed in research and need to be made available in an accessible way with flexible integration into existent workflows.</p>
      </td>
      <td>Geoscience, Geoinformatics, Geodesy</td>
      <td>Höfle group</td>
      <td>
        <p></p>
      </td>
      <td><a href="https://github.com/3dgeo-heidelberg/py4dgeo">https://github.com/3dgeo-heidelberg/py4dgeo</a></td>
    </tr>
    <tr>
      <td>
        <p><strong>pynbody</strong></p>
        <p>pynbody is an analysis package for astrophysical N-body and Smooth Particle Hydrodynamics simulations.</p>
      </td>
      <td>Physics, Astrophysics</td>
      <td>Buck group</td>
      <td>
        <p>ascl code record ascl:1305.002</p>
        <p>10.1093/mnras/stt788</p>
      </td>
      <td><a href="https://pynbody.github.io/">https://pynbody.github.io/</a></td>
    </tr>
    <tr>
      <td>
        <p><strong>pytreedb</strong></p>
        <p>pytreedb is a Python software package providing an object-based library to provide a simple database interface and REST API of vegetation tree objects that were captured as 3D point clouds. The main objective is to provide a Python library for the storage and sharing of single tree-based point clouds and all relevant (forest inventory) tree measurements. The tree data include all tree related information, measurements, metadata, geoinformation and also links to the 3D point clouds linked in any file format (e.g. LAS/LAZ). GeoJSON including all tree-related information is used as data format for data exchange and visualization with other software (e.g. direct import into most GIS). Thereby view and modification of the tree datasets (*.geojson) is straightforward.</p>
      </td>
      <td>Geoscience, Geoinformatics, Geodesy</td>
      <td>Höfle group</td>
      <td>
        <p>10.5281/zenodo.7551310</p>
      </td>
      <td><a href="https://github.com/3dgeo-heidelberg/pytreedb">https://github.com/3dgeo-heidelberg/pytreedb</a></td>
    </tr>
    <tr>
      <td>
        <p><strong>RCarb</strong></p>
        <p>RCarb is a package for dose rate modelling of carbonate-rich samples in the context of trapped charged dating (e.g., luminescence dating) applications.</p>
      </td>
      <td>Geoscience, Geoinformatics, Applied Physics</td>
      <td>Kreutzer group</td>
      <td>
        <p>Ancient TL, 37(2), 1-8 (2019)</p>
      </td>
      <td><a href="https://cran.r-project.org/web/packages/RCarb/index.html">https://cran.r-project.org/web/packages/RCarb/index.html</a></td>
    </tr>
    <tr>
      <td>
        <p><strong>RADMC-3D</strong></p>
        <p>RADMC-3D is a code package for diagnostic radiative transfer calculations in astronomy and astrophysics. It calculates, for a given geometrical distribution of gas and/or dust, what its images and/or spectra look like when viewed from a certain angle, allowing modelers to compare their models with observed data.</p>
      </td>
      <td>Physics, Astronomy, Astrophysics</td>
      <td>Dullemond group</td>
      <td>
        <p>ads bibcode 2012ascl.soft02015D</p>
      </td>
      <td><a href="https://www.ita.uni-heidelberg.de/~dullemond/software/radmc-3d/index.php">https://www.ita.uni-heidelberg.de/~dullemond/software/radmc-3d/index.php</a></td>
    </tr>
    <tr>
      <td>
        <p><strong>RAMD</strong></p>
        <p>The Random Acceleration Molecular Dynamics (RAMD) method can be used to carry out molecular dynamics simulations with an additional randomly oriented force applied to a molecule in the system. Originally this was implemented from the MCM group in Amber 8 (not maintained). Recently the Amber group has provided the functionality integrated in Amber 20. Here is an implementation of the method in NAMD and GROMACS.</p>
      </td>
      <td>Biology, Bioinformatics,</td>
      <td>Wade group</td>
      <td>
        <p>10.1006/jmbi.2000.4154</p>
        <p>10.1007/s008940050053</p>
        <p>10.1021/acs.jctc.8b00230</p>
        <p>10.3389/fmolb.2019.00036</p>
      </td>
      <td><a href="https://www.h-its.org/downloads/ramd/">https://www.h-its.org/downloads/ramd/</a></td>
    </tr>
    <tr>
      <td>
        <p><strong>RASPD+</strong></p>
        <p>RASPD+ (RApid Screening of hit molecules for target proteins via Physicochemical Descriptors+) is a computationally fast protocol for identifying lead-like molecules based on predicted binding freeenergy against a target protein with a 3D structure and a defined ligand binding pocket.&nbsp;</p>
      </td>
      <td>Biology, Bioinformatics, Machine Learning</td>
      <td>Wade group</td>
      <td>
        <p>https://doi.org/10.26434/chemrxiv.12636704.v2</p>
      </td>
      <td><a href="https://github.com/HITS-MCM/RASPDplus">https://github.com/HITS-MCM/RASPDplus</a></td>
    </tr>
    <tr>
      <td>
        <p><strong>RLumCarlo</strong></p>
        <p>RLumCarlo is a collection of functions to simulate luminescence production in dosimetric materials using Monte Carlo methods. Implemented are models for delocalised transitions, localised transitions and tunnelling transitions. Supported stimulation methods are thermal luminescence (TL), continuous-wave optically stimulated luminescence (CW-OSL), linearly-modulated optically stimulated luminescence (LM-OSL), linearly-modulated infrared stimulated luminescence (LM-IRSL), and isothermal luminescence (ITL or ISO-TL).</p>
      </td>
      <td>Geoscience, Geoinformatics, Applied Physics</td>
      <td>Kreutzer group</td>
      <td>
        <p>10.32614/RJ-2021-043</p>
      </td>
      <td><a href="https://cran.r-project.org/web/packages/RLumCarlo/index.html">https://cran.r-project.org/web/packages/RLumCarlo/index.html</a></td>
    </tr>
    <tr>
      <td>
        <p><strong>rxylib</strong></p>
        <p>rxylib provides access to the 'xylib' C++ library to import xy data from powder diffraction, spectroscopy and other experimental methods.</p>
      </td>
      <td>Geoscience, Geoinformatics, Applied Physics</td>
      <td>Kreutzer group</td>
      <td>
        <p></p>
      </td>
      <td><a href="https://cran.r-project.org/web/packages/rxylib/index.html">https://cran.r-project.org/web/packages/rxylib/index.html</a></td>
    </tr>
    <tr>
      <td>
        <p><strong>SCOOP template engine</strong></p>
        <p>The purpose of the SCOOP template engine is to facilitate the preparation of manuscripts in LaTeX for publication in scientific journals. It allows the user to concentrate on the content, rather than the layout. The layout, which depends on the journal, will be automatically generated.&nbsp;</p>
      </td>
      <td>&nbsp;</td>
      <td>Herzog group</td>
      <td>
        <p></p>
      </td>
      <td><a href="https://pypi.org/project/scoop-template-engine/">https://pypi.org/project/scoop-template-engine/</a></td>
    </tr>
    <tr>
      <td>
        <p><strong>SDA / webSDA</strong></p>
        <p>SDA7 can be used to carry out Brownian dynamics simulations of the diffusional association in a continuum aqueous solvent of two solute molecules, e.g. proteins, or of a solute molecule to an inorganic surface. SDA7 can also be used to simulate the diffusion of multiple proteins, in dilute or concentrated solutions, e.g., to study the effects of macromolecular crowding.</p>
      </td>
      <td>Biology, Bioinformatics</td>
      <td>Wade group</td>
      <td>
        <p>10.1002/jcc.23971</p>
        <p>10.1006/meth.1998.0588</p>
        <p>10.1016/S0006-3495(97)78838-6</p>
      </td>
      <td><a href="https://mcm.h-its.org/sda/doc/doc_sda7/index.html">https://mcm.h-its.org/sda/doc/doc_sda7/index.html</a></td>
    </tr>
    <tr>
      <td>
        <p><strong>SEVN</strong></p>
        <p>SEVN (Stellar EVolution for 𝑁-body) is a rapid binary population synthesis code. It takes as input the initial conditions of stars or binaries (masses, spin, semi-major axis, eccentricity) and evolves them. SEVN calculates stellar evolution by interpolating pre-computed sets of stellar tracks. Binary evolution is implemented by means of analytic and semi-analytic prescriptions.</p>
      </td>
      <td>Physics, Astrophysics</td>
      <td>Mapelli group</td>
      <td>
        <p>10.1093/mnras/stad1630</p>
        <p>10.3847/1538-4357/ab584d</p>
      </td>
      <td><a href="https://gitlab.com/sevncodes/sevn">https://gitlab.com/sevncodes/sevn</a></td>
    </tr>
    <tr>
      <td>
        <p><strong>Sketch Map Tool</strong></p>
        <p>The Sketch&nbsp;Map Tool is an easy-to-use tool for participatory sketch mapping through offline collection, digitization and georeferencing of local spatial knowledge. The tool has a variety of applications. For example, do you want to work together with people in a community to map their experience and perception of risk in their neighbourhood in a paper-based format, but still be able to quickly analyse the results digitally? Then, the Sketch&nbsp;Map Tool is exactly what you need!</p>
      </td>
      <td>Geoscience, Geoinformatics</td>
      <td>Zipf group</td>
      <td>
        <p></p>
      </td>
      <td><a href="https://github.com/GIScience/sketch-map-tool">https://github.com/GIScience/sketch-map-tool</a></td>
    </tr>
    <tr>
      <td>
        <p><strong>Spatial Model Editor</strong></p>
        <p>A user friendly GUI editor to create and edit spatial SBML models of bio-chemical reactions and simulate them using the dune-copasi solver for reaction-diffusion systems.</p>
      </td>
      <td>Biology, Bioinformatics</td>
      <td>Kummer group</td>
      <td>
        <p>10.5281/zenodo.8316368</p>
      </td>
      <td><a href="https://spatial-model-editor.github.io/">https://spatial-model-editor.github.io/</a></td>
    </tr>
    <tr>
      <td>
        <p><strong>SYCAMORE</strong></p>
        <p>SYCAMORE is a system that provides you with a facilitated access to a number of tools and methods in order to build models of biochemical systems, view, analyse and refine them, as well as perform quick simulations. SYCAMORE is not intended to substitute for expert simulation and modeling software packages, but might interact with those. It is rather intended to support and guide system biologists when doing computational research. One important function of SYCAMORE is to allow you to build a draft model of your system of interest in such a way that kinetic expressions and parameters are as close to reality as possible. We want to emphasize that the resulting model still has a draft character and should not be taken as "the final model". However, setting up your model in such a way that parameters etc. are as close to reality as possible on the basis of literature data and computational parameter estimation methods should faciliate any parameter fitting methods that you want to employ later on.</p>
      </td>
      <td>Biology, Bioinformatics</td>
      <td>Wade group</td>
      <td>
        <p>10.1093/bioinformatics/btn207</p>
      </td>
      <td><a href="http://sycamore.h-its.org/sycamore/">http://sycamore.h-its.org/sycamore/</a></td>
    </tr>
    <tr>
      <td>
        <p><strong>TiSR</strong></p>
        <p>To improve upon and accelerate the thermodynamic equation of state development process, we introduce thermodynamics-informed symbolic regression (TiSR). It aims to combine a symbolic regression base with the extensions required to work with often strongly scattered experimental data, different residual pre- and post-processing options, and additional features required to consider thermodynamic equations of state development. The project and the code are a work-in-progress and there are many more features planned.</p>
      </td>
      <td>Machine learning in thermodynamics, explainable AI</td>
      <td>Herzog group/Viktor Martinek</td>
      <td>
        <p>10.48550/arXiv.2309.02805, 10.5281/zenodo.8317547</p>
      </td>
      <td><a href="https://github.com/scoop-group/TiSR">https://github.com/scoop-group/TiSR</a></td>
    </tr>
    <tr>
      <td>
        <p><strong>TRAPP</strong></p>
        <p>TRAnsient Pockets in Proteins with Druggability Calculation allows analysis of transient binding pockets and druggability indices in proteins. The TRAPP webserver, and the downloadable command line version, is intended to aid the discovery of ligands that bind in transient subpockets in proteins. The TRAPP webserver offers a range of tools to explore binding site motions ranging from local side chain fluctuations to global backbone motions; analysis of the binding site dynamics in simulated protein motion trajectories or trajectories provided by the user; a tool for tracking, analysis and visualization of protein cavity dynamics in an ensemble of protein structures or in protein trajectories; detection of transient pockets and subpockets that may appear due to protein motion; combined analysis of the protein sequence conservation and the variations of the physicochemical properties of a binding pocket; assessment of the binding pocket druggability and its variation in simulated trajectories and visualization of the pocket regions that favorably contribute to the pocket druggability. TRAPP is not designed to identify all of a protein's binding pockets, but rather to trace changes in the spatial and physicochemical properties of a specified pocket in a protein that may arise due to the protein's flexibility.</p>
      </td>
      <td>Biology, Bioinformatics</td>
      <td>Wade group</td>
      <td>
        <p>10.1021/acs.jcim.9b01185</p>
      </td>
      <td><a href="https://trapp.h-its.org/">https://trapp.h-its.org/</a></td>
    </tr>
    <tr>
      <td>
        <p><strong>UNICORN</strong></p>
        <p>UNICORN uses a spatially resolved simulation area in which an Active Galactic Nuclei jet is divided into disks along its axis. In each disc, a photon population and a forward and a backward electron and proton population are considered. These move along the jet axis by advection. In parallel, a given velocity profile of the background plasma is considered. After a certain time, due to the scattering of this plasma, an equilibrium is established between the particles moving forward and backward. For a shock in the background plasma, there is then a direct Fermi-I acceleration in the model.</p>
      </td>
      <td>Physics, Astrophysics</td>
      <td>Spanier group</td>
      <td>
        <p>10.1051/0004-6361/201424159</p>
        <p>10.1016/j.astropartphys.2018.02.008</p>
      </td>
      <td><a href="https://www.ita.uni-heidelberg.de/~fspanier/software/unicorn.shtml?lang…">https://www.ita.uni-heidelberg.de/~fspanier/software/unicorn.shtml?lang…</a></td>
    </tr>
    <tr>
      <td>
        <p><strong>UTOPIA</strong></p>
        <p>Utopia is a comprehensive modelling framework for complex and adaptive systems that provides an integrated, automated workflow from data generation to data analysis and visualization. It thereby uses C++ for model implementation and Python for data analysis and is designed for user-friendliness and to facilitate collaborative research.</p>
      </td>
      <td>Physics, Complexity Science</td>
      <td>Roth group</td>
      <td>
        <p>10.21105/joss.02165</p>
        <p>10.1007/978-3-030-50436-6_32</p>
      </td>
      <td><a href="https://utopia-project.org">https://utopia-project.org</a></td>
    </tr>
    <tr>
      <td>
        <p><strong>VOSTOK</strong></p>
        <p>VOSTOK (Voxel Octree Solar Toolkit) is a command-line tool to compute a detailed model of incoming solar radiation distribution on a patch of land, including structures like buildings and vegetation, represented by a 3D point cloud data set. "Vostok" is also the Russian word for "east" - the direction in which the sun rises.</p>
      </td>
      <td>Geoscience, Geoinformatics, Geodesy</td>
      <td>Höfle group</td>
      <td>
        <p>10.11588/data/QNA02B</p>
      </td>
      <td><a href="https://github.com/3dgeo-heidelberg/vostok">https://github.com/3dgeo-heidelberg/vostok</a></td>
    </tr>
    <tr>
      <td>
        <p><strong>VSTT</strong></p>
        <p>Visuomotor Serial Targetting Task (VSTT) is an open source Python GUI tool for designing, running and analyzing motor skill acquisition experiments.</p>
      </td>
      <td>Sports Science</td>
      <td>Wanner group/SSC</td>
      <td>
        <p></p>
      </td>
      <td><a href="https://github.com/ssciwr/vstt">https://github.com/ssciwr/vstt</a></td>
    </tr>
    <tr>
      <td>
        <p><strong>XLUM</strong></p>
        <p>XLUM is an open data format for exchange and long-term preservation of luminescence data.</p>
      </td>
      <td>Geoscience, Geoinformatics, Applied Physics</td>
      <td>Kreutzer group</td>
      <td>
        <p>10.5194/gchron-5-271-2023</p>
      </td>
      <td><a href="https://doi.org/10.5281/zenodo.7362363">https://doi.org/10.5281/zenodo.7362363</a></td>
    </tr>
  </tbody>
</table>
```
