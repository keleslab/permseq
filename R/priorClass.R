setClassUnion("OptionalNumeric",c("numeric","NULL"))
setClassUnion("OptionalCharacter",c("character","NULL"))
setClassUnion("OptionalList",c("list","NULL"))
setClassUnion("OptionalCharacterNumeric",c("character","numeric","NULL"))
setClass("Prior",representation=representation(
                   dnaseName = "OptionalCharacterNumeric",
                   dnaseAlign = "OptionalList",
                   dnaseKnots = "OptionalNumeric",
                   dnaseThres = "OptionalNumeric",
                   posLoc_bychr = "OptionalList",
                   dnaseHistone = "OptionalList",

                   histoneName = "OptionalCharacter",
                   histoneNum = "OptionalNumeric",
                   histoneAlign = "OptionalList",
                   histoneGrpL = "OptionalCharacter",

                   chipName = "OptionalCharacter",
                   chipNum = "OptionalNumeric",
                   chipAlign = "OptionalList",
                   chipSAM = "OptionalCharacter",
                   chipAllocate = "OptionalCharacter",
                   chipUni = "OptionalCharacter",
                   chipAlignFormat = "OptionalCharacter",

                   dataNum = "OptionalNumeric",
                   fragL = "OptionalNumeric",
                   chrList = "OptionalCharacter",
                   bwaInfo = "OptionalList",
                   bowtieInfo = "OptionalList",
                   csemDir = "OptionalCharacter",
                   picardDir = "OptionalCharacter",
                   outfileLoc = "OptionalCharacter",
                   prior = "OptionalCharacter",
                   chrom.ref = "OptionalCharacter"
 
  )
)

















