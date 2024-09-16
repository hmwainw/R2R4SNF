# Reprocessing Data README

In fuel cycles were used fuel is reprocessed, the composition of the waste streams must be known to characterize the management and disposal of those wastes. The composition of the waste streams is determined by the type of used fuel, the selected process, the choice of recovered product, and the separation efficiency. In this README, descriptions for each reprocessing method are provided and, where applicable, justification is given for the assumed partition of elements among waste streams. The reprocessing data is contained in the .sep data files.

| Data file                            | Wastes       | Products     | 
|--------------------------------------|:------------:|:------------:|
| `echem_th-u_fces.sep`                | $\checkmark$ |              |
| `echem_th-u-pu_fces.sep`             | $\checkmark$ |              |
| `echem_th-u-tru_fces.sep`            | $\checkmark$ |              |
| `echem_u-pu_fces.sep`                | $\checkmark$ | $\checkmark$ |
| `echem_u-tru_fces.sep`               | $\checkmark$ | $\checkmark$ |
| `meltrefining_th-u-tru-fp_fces.sep ` | $\checkmark$ |              |
| `meltrefining_u-tru-fp_fces.sep `    | $\checkmark$ |              |
| `thorex_th-u-tru_fces.sep`           | $\checkmark$ |              |
| `thorex_th-u_fces.sep`               | $\checkmark$ |              |
| `urex_th-u-tru_fces.sep`             | $\checkmark$ |              |
| `urex_th-u_fces.sep`                 | $\checkmark$ |              |
| `urex_u-pu_cs-sr.sep`                | $\checkmark$ | $\checkmark$ |
| `urex_u-pu_fces.sep`                 | $\checkmark$ | $\checkmark$ |
| `urex_u-tru_cs-sr.sep`               | $\checkmark$ | $\checkmark$ |
| `urex_u-tru_fces.sep`                | $\checkmark$ | $\checkmark$ |
| `urex_u_cs-sr.sep`                   | $\checkmark$ | $\checkmark$ |
| `urex_u_fces.sep`                    | $\checkmark$ | $\checkmark$ |


## Electrochemical reprocessing

Metal fuel: produces metal and ceramic wastes, all volatile FP as gas (“echem”)

Oxide fuel: produces metal and ceramic wastes, all volatile and some semi-volatile FP in gas (“pyrox”)

### FCE&S Description: 

Electro-chemical processing is a compact efficient method developed for recycling fast reactor metallic fuel on-site with very short cooling times. This process has been demonstrated during the U.S. fast reactor program and is being applied at the Fuel Conditioning Facility in INL for treating EBR-II used fuel. The method provides separation of actinide elements from fission products by means of an electrorefining step. The process of electrorefining is based on well-understood electrochemical concepts. A similar technology, requiring an additional head-end processing step, can be applied with equal success to fuel types other than metallic fuel, such as light-water reactor oxide fuel (PYROX).

### Partitioning:

The elements in the used fuel can be partitioned into four groups. The first group is gaseous fission products, such as the noble gases and carbon-14. These are collected during chopping and electrorefining in the off gas system and stored in canisters. Of the remaining three groups that do not form gases, two are relatively easy to classify into one of the two waste streams based on the free energy of formation of their chlorides.

Rather than reacting with the salt, noble-metal fission products remain with cladding hulls in the anode basket. The free energies of formation for chlorides of these elements are shown on the right hand columns of the table below. Reactive fission products, such as alkali-metal, alkaline-earth, some rare earth, and halide fission products form chloride compounds in the salt. These elements are shown in the first two columns on the left in the table below. The third group includes any remaining rare earth (typically lanthanide) fission products, Zr, and all actinides. These elements exist in equilibrium as metals (possibly dissolved in cadmium) or as chlorides in the salt. These elements are shown in the middle two columns of the table.

#####Table 1. Free energies of formation of SNF chlorides $\Delta G_f^0(500 ^oC)$ in units of kcal/mole [1-3]

|Element  |$\Delta G_f^0$|Element  |$\Delta G_f^0$|Element  |$\Delta G_f^0$| 
|---------|--------------|---------|--------------|---------|--------------|
|MnCl$_2$ | -90.6        |CmCl$_3$ | -64.0        |TlCl     | -37.4        |
|BaCl$_2$ | -87.9        |PuCl$_3$ | -62.4        |InCl     | -34.9        |
|CsCl     | -87.8        |AmCl$_3$ | -62.1        |CdCl$_2$ | -32.3        |
|RbCl     | -87.0        |NpCl$_3$ | -58.1        |FeCl$_2$ | -29.9        |
|KCl      | -86.7        |PbCl$_2$ | -58.1        |CuCl     | -27.9        |
|SrCl$_2$ | -84.7        |UCl$_3$  | -55.2        |NbCl$_5$ | -26.7        |
|LiCl     | -82.5        |CoCl$_2$ | -49.3        |CuCl$_2$ | -25.7        |
|NaCl     | -81.2        |ZrCl$_4$ | -46.6        |MoCl$_4$ | -16.8        |
|CaCl$_2$ | -80.7        |NiCl$_2$ | -45.1        |TcCl$_4$ | -11.0        |
|ZnCl$_2$ | -74.7        |         |              |RhCl$_3$ | -10.0        |
|CrCl$_2$ | -71.2        |         |              |PdCl$_2$ | -9.0         |
|LaCl$_2$ | -70.2        |         |              |RuCl$_4$ | -6.0         |
|PrCl$_2$ | -69.0        |         |              |         |              |
|CeCl$_2$ | -68.6        |         |              |         |              |
|NdCl$_2$ | -67.9        |         |              |         |              |
|YCl$_3$  | -65.1        |         |              |         |              |


Based on their stability in the salt [1-3], Table 2 identifies fission product elements that are assumed to be partitioned 100% into the salt and metal waste streams. This partition excludes the rare earth elements with less stable chlorides and the actinides (will be dealt with later). Normally, due to the similarity of their chemical behavior with that of the actinides, some amount of rare earths will typically be transported along with actinides to the liquid cadmium cathode, where uranium and transuranics accumulate. In this analysis, because the separation efficiency is predefined from the FCE&S study, it is assumed that no rare earth fission products are found in any of the product streams. Therefore, in this analysis, they are left to distribute between the salt and liquid cadmium in the electrorefiner.

#####Table 2. Partitioning of elements in electrorefiner wastes

| Salt (ceramic waste form) | Metal (metal waste form) |
|---------------------------|--------------------------|
| Alkali metals (Li, Na, K, Rb, Cs, Fr) | Noble metal FP (Y, Nb, Mo, Tc, Ru, Rh, Pd, Ag, Cd, Sn, Te, Zr, etc.) |
| Alkaline metal FP and decay products (Be, Mg, Ca, Sr, Ba, Ra) | Cladding materials (Fe, Cr, Ni, Zr) |
| Lanthanides with stable chlorides (La, Ce, Pr, Nd, Sm, Eu) | FP with anions that will dissolve in the salt (Br, I, Se, Te, As, Sb) |

In the FCE&S study, it is assumed that unrecovered actinides are sent for disposal. For some fuel cycles, only U and Pu are recycled, while the remaining TRU is treated as waste. However, in the pyroprocess, all transuranics are recovered to some extent in the liquid cadmium cathode. Therefore, the selective recovery of actinides by pyroprocessing (say, separating Pu from Am) may require extra steps or special operation of the electrorefiner. The determination of electrorefiner operation is outside the scope of this work. For this analysis, it is assumed that all actinides not recovered as product remain in the salt phase [4]. Additionally, the remaining rare earth elements are assumed to remain in the salt as well; both should end up in the ceramic waste form. 

The justification for this assumption is as follows: the composition of the U/TRU (or U/Pu) product mixture extracted at the liquid cathode depends on the amounts of U, Pu, rare earth metals, and cadmium in the electrorefiner anode, salt, and cathode. There exists a salt composition that can result in recovered product of any composition. The composition of the salt is controlled by the addition or removal of species during electrorefiner operation [5].

For pyroprocessing of oxide fuels, the same justification is applied for the partitioning of elements between the salt and metal waste streams. However, due to the additional head end processing, some semi-volatile fission products will be released into the off gas collection system. Thus, it is assumed (rather arbitrarily) that 20% of the semi-volatile FP (Se, Br, Rb, Mo, Tc, Cd, Sb, I, Cs) end up in the off gas waste along with the fully volatile FP. The remaining semi-volatile FP are partitioned as specified above.

###### References
1. National Research Council. “Electrometallurgical Techniques for DOE Spent Fuel Treatment – Final Result”. National Academy Press, Washington, D.C. (2000).

2. L. Pancratz. “Thermodynamic Properties of Halides”, United States Department of the Interior, Bureau of Mines, Bulletin 674.

3. Fuger, J., V. Parker, W. Hubbard and F. Oetting. “The Chemical Thermodynamics of Actinide Elements and Compounds, Part 8, The Actinide Halides”, International Atomic Energy Agency, Vienna, Austria (1983).

4. K. M. GOFF, J. C. WASS, K. C. MARSDEN, and G. M. TESKE. “Electrochemical Processing of Used Nuclear Fuel”. Nuclear Engineering and Technology, v43 n4, pp. 335-342 (2011).

5. J. P. Ackerman. “Chemical Basis for Pyrochemical Reprocessing of Nuclear Fuel”. Industrial & Engineering Chemistry Research, v30 n1, pp. 141-145 (1991)

6. S. Frank, W. Ebert, B. Riley, H. S. Park, Y. Z. Cho, C. H. Lee, M. K. Jeon, J. H. Yang, H. C. Eun. "Waste Stream Treatment and Waste Form Fabrication for Pyroprocessing of Used Nuclear Fuel" INL/EXT-14-34014 (2015).


## UREX repcrocessing

UREX+1: Produces U, TRU, Ln-FP, non-Ln-FP

UREX+2: Produces U, Pu/Np, Am/Cm/Ln, non-Ln FP

UREX+3: Produces U, Pu/Np, Am/Cm, Ln-FP, non-Ln FP

UREX+4: Produces U, Pu/Np, Am, Cm, Ln-FP, non-Ln FP

### FCE&S Description: 

The generic name UREX+ reflects a design that initially extracts uranium (U) and technetium (Tc), which are subsequently separated by ion exchange, and followed by additional processing to recover transuranic (TRU) and fission product (FP) elements in different combinations. UREX is generally followed by sequestration of cesium and strontium (Cs/Sr). The major variants are then distinguished by the ensuing separation of the balance of the FPs from the TRU elements.

UREX+ was intended primarily to provide options for the processing of spent nuclear fuel without separating a pure plutonium stream. UREX+1 maintains the TRU elements as a group, while the lanthanide (Ln) and non-lanthanide FPs are recovered separately. UREX+2 arranges these into three products: Pu/Np, Am/Cm/Ln, and the remaining non-Ln FPs. UREX+3 also recovers Pu/Np, but further separates Am and Cm from the Ln FPs. UREX+4 adds the separation of Am from Cm to the UREX+3 process.

There are minor variations of each UREX+ technology both in terms of processing and in terms of products generated, specifically variants that recover U with the Pu/Np-bearing streams. All of the UREX+ process variants were intended to provide the products needed for specific applications. For example, the UREX+3c process is intended for LWR recycle fuel, while the UREX+1a process is intended for FR recycle fuel; these two variants have been demonstrated at laboratory scale. Aqueous processes have exhibited very high recovery efficiencies and very high decontamination levels. Though the more complex extraction chemistries of UREX+ have been shown to be feasible, more work is required to develop, demonstrate, and scale these up, along with the associated peripheral processes.

### Partitioning

The UREX+ suite of aqueous separations processes all begin with the UREX (**U**ranium **E**xtraction) process to separate uranium from the rest of spent fuel. Numerous other waste streams are produced, depending in part on the types of separations performed after the initial removal of uranium. UREX+1 and UREX+1a were designed for homogeneous recycle of all transuranics in fast spectrum reactors. UREX+2, +3, and +4 were designed for selective recycling. These additional separations do not significantly improve the repository performance, but affect the ease with which recovered material can be handled. If fuel fabrication requires remote handling, costs will certainly be higher, and the technology required for remote handling is immature [1]. 

In UREX+1 and +1a, the UREX process is used to separate uranium from all other elements in the used fuel. In UREX+2, +3, and +4, uranium is coextracted with the Pu/Np product. The remaining elements go on to additional processes for further separations. In UREX+2, the TRUEX process is used to separate Am+Cm+Ln from the remaining FP. In UREX+3, the Am+Cm+Ln product from TRUEX is treated in the TALSPEAK process, which separates Am+Cm from the lanthanides. Finally, In UREX+4, the Am/Cm stream is separated into Am and Cm streams.

##### Table 1. Product streams from different UREX+ processes [1]

| Product | 1        | 2  | 3 | 4      | 5              | 6        | 7  |
|---------|----------|----|---|--------|----------------|----------|----|
| UREX+1  | U        | Tc | I | Cs, Sr | Other FPs      | TRU+Ln   |    |
| UREX+1a | U        | Tc | I | Cs, Sr | FPs (incl. Ln) | all TRU  |    |
| UREX+2  | U, Pu+Np | Tc | I | Cs, Sr | Other FPs      | Am+Cm+Ln |    |
| UREX+3  | U, Pu+Np | Tc | I | Cs, Sr | FPs (incl. Ln) | Am+Cm    |    |
| UREX+4  | U, Pu+Np | Tc | I | Cs, Sr | FPs (incl. Ln) | Am       | Cm |

Not shown in Table 1 are the fission product gases Kr, Xe, and C-14, which are captured during chopping and dissolution. 

In Table 1, the first four product streams are common to each UREX+ process. The uranium is recovered for reuse or disposal as LLW. Technetium remains undissolved and is recovered along with cladding hulls. Iodine is separated from the off gas system and converted to KI for HLW disposal. Both Tc and I are long-lived FPs that are mobile in most geologic environments and are therefore key contributors to long-term repository dose. The Cs-Sr waste is highly purified and stored to reduce heat emission in repository wastes. 

Where separations are assumed incomplete (for example, if it assumed that 99% of U and Pu is recovered) the remaining material is assumed to end up in a waste stream depending on which materials are desired for recovery. 

* If U is recovered alone, the UREX+1 process will be used and unrecovered U will end up in the TRU/Ln stream (TRUEX raffinate). File: urex\_u\_fces.sep

* In some cases, the FCE&S study specifies thorium recovery as a part of the UREX+ processes. This capability is largely unreported in literature. Therefore, for the purposes of this analysis, if Th is recovered along with uranium, it is assumed to take place during the main UREX step as if the UREX+1 processed was used. File: urex\_u-th\_fces.sep

* If U and TRU are recovered, the UREX+1a, +3, or +4 processes may be used. Unrecovered actinides are disposed of in the waste stream that contains that lanthanides. File: urex\_u-tru\_fces.sep

* If U and Pu are recovered without other TRU, the UREX+2 process will be used, and the unrecovered actinides will be disposed of in the TRU (Am+Cm+Ln) waste stream. File: urex\_u-pu\_fces.sep

<!--Also need to include Th/TRU and Th/TrTh and others?! Maybe do this on a case-by-case basis, have Madison figure this out.-->

###### References
1. J. J. Laidler. "GNEP Spent Fuel Processing; Waste Streams and Disposition Options". Nuclear Waste Technical Review Board, Washington, D.C., 15 May 2007, <http://www.nwtrb.gov/meetings/2007/may/laidler.pdf>.