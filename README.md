# seasonality
This project is for the calculation of seasonality from daily CHIRPS rainfall and Hoffman ETo data

## Outputs

Files are saved in a folder with a name derived from the input parameters, for example:   `../S2mm0-Pad3x3-D1mm25-D2mm20-D2l2-AIt0.5-SqMT-SqMLTT-MxGp1-MiSL1-MaSL1-ClAIF-S2Pr0.25-S2l1-AIsF-RBT-S1AITRUE`  

Where:
  `S2mm` = The minimum rainfall required for second rainy season at a location (`MinRain`)
  `-Pad` = PadBack,"x",PadForward,
  "-D1mm",D1.mm,
  "-D2mm",D2.mm,
  "-D2l",D2.len,
  "-AIt",AI.t, 
  "-SqM",substr(Do.SeqMerge,1,1),
  "-SqMLT",substr(Do.SeqMerge.LT,1,1),
  "-MxGp",MaxGap,
  "-MiSL",MinStartLen,
  "-MaSL",MaxStartSep,
  "-ClAI",substr(ClipAI,1,1),
  "-S2Pr",Season2.Prop,
  "-S2l",MinLength,
  "-AIs",substr(AI_Seasonal,1,1),
  "-RB",substr(RollBack,1,1),
  "-ST",SOSSimThresh*100,
  "-S1AI",S1.AI
