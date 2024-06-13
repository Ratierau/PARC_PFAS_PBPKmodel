# PBK models for PFAS 
This repository contains R code for developing PBK model for describing PFAS biodistribution in multiple species and deploying as web services on Jaqpot (https://www.jaqpot.org/). This work is part of the SCENARIOS EU H2020 project (https://scenarios-project.eu/). The list of available PBK models will be updated throughout the lifespan of the project.

## List of PBK models
- Humans:
    * PFOS/PFOA (Loccisano et al., 2011)
    * PFAAs (Fabrega et al., 2015)
    * PFOS/PFOA (Fabrega et al., 2016)
- Rats:
    * PFOA [male, female] (Loccisano et al., 2012; Worley and Fisher, 2015)
    * PFOS [male, female] (Loccisano et al., 2012)
    * PFHxS [male, female] (Kim et al., 2018)
    * PFNA [male, female] (Kim et al., 2019)
    * PFDA [male, female] (Kim et al., 2019)
- Fish:
   * PFOS Rainbow trout PBPK (Vidal et al., 2020)
   * PFAAs Protein Binding Fish PBK (Ng and Konrad Hungerbühler, 2013)

## Folder Contents
*  Bernstein et al.2021 --> This folder transfers the PFOA, PFOS, PFHxS, PFNA and PFDA PBK code presented in Bernstein et al. (2021) to R code in an appropriate format, so that it can be readily exposed as a web service on Jaqpot.
*  Fabrega et al.2015 --> The file "PFAS_PBK_Fabrega_2015_Jaqpot_upload.R" consists of the model and can be used to upload the mon on Jaqpot as a web service.
*  Fabrega et al.2016 --> The file "PFOS_PFOA Human PBK.R" consists of the model and can be used to upload the mon on Jaqpot as a web service.
*  Ng and Hungerbühler 2013 --> The file "Ng_Hungerbuhler_2013.R" consists of the model and can be used to upload the mon on Jaqpot as a web service.
*  PFAS Rainbow trout PBK --> This folder contains the code of an under development PBK model. The purpose of the model is to predict the concentration of PFAAs substances in Rainbow trout after being dietary exposed.
*  Vidal et al.2019 --> The file "PFOS_Rainbow_trout_PBTK.R" contains the model written in R and the code to easily upload the model on Jaqpot.
*  Vidal et al.2020 --> The file "Jaqpot_final.R" contains the model written in R and the code to easily upload the model on Jaqpot.
*  Worley and Fisher.2015 --> The file "Worley_rat.R" contains the rat model written in R and the code to easily upload the model on Jaqpot.

## Necessary R packages
* `deSolve`  --> A library containing functions for solving ODE systems
* `jaqpotr` --> A library containing functions for uploading models on the Jaqpot server. Installation instructions can be found in https://www.jaqpot.org/docs/r

## References
* Bernstein, A.S., Kapraun, D.F. and Schlosser, P.M. (2021). A Model Template Approach for Rapid Evaluation and Application of Physiologically Based Pharmacokinetic Models for Use in Human Health Risk Assessments: A Case Study on Per- and Polyfluoroalkyl Substances. Toxicological Sciences, 182(2), pp.215–228. doi:https://doi.org/10.1093/toxsci/kfab063.Francesc Fàbrega, Kumar, V., Benfenati, E., Schuhmacher, M., Domingo, J.L. and Nadal, M. (2015). Physiologically based pharmacokinetic modeling of perfluoroalkyl substances in the human body. 97(6), pp.814–827. doi:https://doi.org/10.1080/02772248.2015.1060976.
*	Fàbrega, F., Kumar, V., Benfenati, E., Schuhmacher, M., Domingo, J.L. and Nadal, M. (2016). Physiologically based pharmacokinetic modeling of perfluoroalkyl substances in the human body. Toxicological & Environmental Chemistry, 97, pp.814-827. doi: https://doi.org/10.1080/02772248.2015.1060976
* Fàbrega, F., Nadal, M., Schuhmacher, M., Domingo, J.L. and Kumar, V. (2016). Influence of the uncertainty in the validation of PBPK models: A case-study for PFOS and PFOA. Regulatory Toxicology and Pharmacology, 77, pp.230–239. doi:https://doi.org/10.1016/j.yrtph.2016.03.009
*	Kim, S.-J., Choi, E.-J., Choi, G.-W., Lee, Y.-B. and Cho, H.-Y. (2018a). Exploring sex differences in human health risk assessment for PFNA and PFDA using a PBPK model. Archives of Toxicology, 93(2), pp.311–330. doi:https://doi.org/10.1007/s00204-018-2365-y
*	Kim, S.-J., Shin, H., Lee, Y.-B. and Cho, H.-Y. (2018b). Sex-specific risk assessment of PFHxS using a physiologically based pharmacokinetic model. Archives of Toxicology, 92(3), pp.1113–1131. doi:https://doi.org/10.1007/s00204-017-2116-5
*	Loccisano, A.E., Campbell, J.L., Andersen, M.E. and Clewell, H.J. (2011). Evaluation and prediction of pharmacokinetics of PFOA and PFOS in the monkey and human using a PBPK model. Regulatory Toxicology and Pharmacology, 59(1), pp.157–175. doi:https://doi.org/10.1016/j.yrtph.2010.12.004
*	Loccisano, A.E., Campbell, J.L., Butenhoff, J.L., Andersen, M.E. and Clewell, H.J. (2012). Comparison and evaluation of pharmacokinetics of PFOA and PFOS in the adult rat using a physiologically based pharmacokinetic model.  Reproductive Toxicology, 33(4), pp.452–467. doi:https://doi.org/10.1016/j.reprotox.2011.04.006
*	Ng, C.A. and Konrad Hungerbühler. (2013). Bioconcentration of Perfluorinated Alkyl Acids: How Important Is Specific Binding? Environmental Science & Technology, 47(13), pp.7214–7223. doi:https://doi.org/10.1021/es400981a
*	Vidal, A., Garric, J., Garric, J. and Rémy Beaudouin. (2020). Temperature effect on perfluorooctane sulfonate toxicokinetics in rainbow trout (Oncorhynchus mykiss): Exploration via a physiologically based toxicokinetic model. Aquatic Toxicology, 225, pp.105545–105545. doi:https://doi.org/10.1016/j.aquatox.2020.105545
*	Worley, R.R. and Fisher, J. (2015). Application of physiologically-based pharmacokinetic modeling to explore the role of kidney transporters in renal reabsorption of perfluorooctanoic acid in the rat. Toxicology and Applied Pharmacology, 289, pp.428–441. doi: https://doi.org/10.1016/j.taap.2015.10.017
