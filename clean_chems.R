library(tidyverse)


uncurated_chemicals <- read_csv("uncurated_chemicals_2023-05-16_12-43-41.csv") %>%
  distinct(raw_cas, raw_chem_name)

#' Extract CASRN numbers from a text string
#'
#' This function takes a text string and extracts all occurrences of a CASRN 
#' number, removing leading zeros.
#'
#' @param x A character vector containing text to be searched for CASRNs
#' @return A list of CASRN numbers found in the text, or NA if none are found
#' @import stringr
#' @import dplyr
casrn_split <- function(x) {
  if (is.character(x)) {
    s <- x %>%
      str_extract_all("[1-9][0-9]{1,6}\\-[0-9]{2}\\-[0-9]", simplify = FALSE)
    
    if (length(s) < 1) {
      s <- NA
    }
  } else {
   s <- NA
  }
  return(s)
}


#' Check that the last digit in the CAS-RN is a valid digit.
#' One way to ensure if a CAS-RN is invalid.
#'
#' @param x A string representing a CAS-RN.
#'
#' @return A boolean: TRUE if checksum is valid, FALSE otherwise.
#'
#' @examples
#' checksum("50-00-0")
#' checksum("123-45-6")

checksum <- function(x) {
  if (!is.character(x)) {
    return(FALSE)
  }
  
  cas <- x %>%
    substr(nchar(x) - 3, 1) %>%
    gsub("-", "", .) %>%
    strsplit(split = "") %>%
    unlist()
  
  q <- 0
  for (i in seq_along(cas)) {
    q <- q + (i) * as.integer(cas[i])
  }
  
  if (q %% 10 == as.integer(substr(x, nchar(x), nchar(x)))) {
    return(TRUE)
  } else {
    return(FALSE)
  }
}

{
  element_pattern <- paste0(
    "(H|He|Li|Be|B|C|N|O|F|Ne|Na|Mg|Al|Si|P|S|Cl|Ar|K|Ca|Sc|Ti|V|Cr|Mn|Fe|Co|Ni|Cu|Zn|Ga|Ge|As|Se|Br|Kr|Rb|Sr|Y|Zr|Nb|Mo|Tc|Ru|Rh|Pd|Ag|Cd|In|Sn|Sb|Te|I|Xe|Cs|Ba|La|Ce|Pr|Nd|Pm|Sm|Eu|Gd|Tb|Dy|Ho|Er|Tm|Yb|Lu|Hf|Ta|W|Re|Os|Ir|Pt|Au|Hg|Tl|Pb|Bi|Po|At|Rn|Fr|Ra|Ac|Th|Pa|U|Np|Pu|Am|Cm|Bk|Cf|Es|Fm|Md|No|Lr|Rf|Db|Sg|Bh|Hs|Mt|Ds|Rg|Cn|Nh|Fl|Mc|Lv|Ts|Og)"
  )
  
  # Define the regex pattern for numbers
  num_pattern <- "(?:[1-9]\\d*)?"
  
  # Define the regex pattern for element groups
  element_group_pattern <- paste0(
    "(?:", element_pattern, num_pattern, ")+"
  )
  
  # Define the regex pattern for parentheses groups
  element_parentheses_group_pattern <- paste0(
    "\\(", element_group_pattern, "\\)", num_pattern
  )
  
  # Define the regex pattern for square bracket groups
  element_square_bracket_group_pattern <- paste0(
    "\\[(?:", element_group_pattern, "|", element_parentheses_group_pattern, ")+\\]", num_pattern
  )
  
  # Combine all patterns into the final regex
  final_regex <- paste0(
    "^(", element_square_bracket_group_pattern, "|", element_parentheses_group_pattern, "|", element_group_pattern, ")+$"
  )
  
  stringr::str_detect("(CH3)2CFCOO(CH2)2Si[NO3(CH3)2]2", final_regex)
}
