################# Splitting the string based on special characters ###########################

paste(unlist(strsplit(x,split = "special_character"))[1:4],collapse = "_")

######## Examples:

paste(unlist(strsplit("dd_smedV4_11074_0_467",split = "_"))[1:4],collapse = "_")
"dd_smedV4_11074_0"

paste(unlist(strsplit("cth1_Trc_comp15661_c0_seq1",split = "_"))[1:4],collapse = "_")
"cth1_Trc_comp15661_c0"