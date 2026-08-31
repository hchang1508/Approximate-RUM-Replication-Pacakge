# Auxiliary functions used by the Fishing-data fit exercises:
# compute.emp.CP, compute.logitCP and dist.prob.
#
# dist.prob is the plain Euclidean distance between choice-probability vectors,
# the metric both fit tables report.

library(dplyr)   # To use combine



# Function to compute empirical choice probabilities
compute.emp.CP <- function(dat){
  # Input:  dat (data.frame, with the variable "alt" being the names of alternatives
  #              and the variable "mode" being the choice indicator)
  # Output: emp.CP (the vector of empirical choice probabilities)
  
  dat.gp <- group_by(dat, alt)    # Group data by alternatives to compute CP
  emp.CP <- summarize(dat.gp, mean = mean(mode)) 
  emp.CP <- emp.CP[, 'mean'] %>% as.matrix() %>% as.numeric()
  
  return(emp.CP)
}





# Function to compute logit predictions given a vector of coefficients and a choice set
compute.logitCP <- function(dat.input, pol.label, b, cs, b.scale = 1){
  # Input:  dat.input (data.frame, with each row corresponding to a distinct alternative.
  #                    It must have a variable alt, indicating the name of alternatives.
  #                    This function takes as an input a dataset which consists of a unique 
  #                    vector of characteristics values for each alternative, such as 
  #                    data on average characteristics.
  #                    When applying this function to data with interactive characteristics,
  #                    apply it to individual- or type-specific characteristics matrix, 
  #                    where type is defined by distinct values of interactive characteristics.)
  #         pol.label (the column names of polynomial terms)
  #         b (numeric vector of coefficients)
  #         cs (choice set, which is a subset of the universal choice set)
  #         b.scale (scale of coefficients)
  # Output: CP (choice probability vector predicted by the vector b*b.scale on the choice set cs)
  # Note: This function computes logit predictions, not mixed-logit ones.
  
  lin.idx    <- as.matrix(dat.input[, pol.label]) %*% b * b.scale  
  lin.idx.cs <- lin.idx[ which( dat.input$alt %in% cs), ]
  CP <- c()
  for (alt in cs){
    lin.idx.diff <- lin.idx.cs - lin.idx[ which(dat.input$alt == alt), ] 
    CP[alt] <- 1/sum(exp(lin.idx.diff))
  }
  
  return(CP)
}


# Function to compute the Euclidean distance between two CP vectors
dist.prob <- function(cp1, cp2) dist( rbind(cp1, cp2) )
# Input:  cp1 (numeric vector of choice probabilities)
#         cp2 (numeric vector of choice probabilities)
# Output: distance between the two CP vectors
