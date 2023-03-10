#' Online Survey on "Exams and Written Papers"
#'
#' The goal of the survey was to estimate the prevalence of various forms of student misconduct such as plagiarizing or cheating in exams.
#' Because students might be reluctant to reveal information on such behaviors, special techniques for sensitive questions were employed in addition to direct questioning.
#' Respondents were randomly assigned to direct questioning or one of five different sensitive question techniques.
#' The dataset contains the (randomized or direct) responses from 4281 students of the University of Bern and ETH Zurich.
#' Each row holds the response to one question for one respondent.
#' The variables are as follows:
#'
#' \itemize{
#'   \item id. Identification code of the respondent
#'   \item RR_response. Binary randomized or direct response
#'   \item Question. Which question was asked
#'   \item expcond. Experimental condition
#'   \item protect. Level of respondent protection
#'   \item subgroup. Subgroups for balanced assignment to experimental conditions
#'   \item sample. Sample group
#'   \item survey duration. Total time to complete survey (in seconds)
#'   \item mobile. Respondent used mobile device (at start of interview)
#'   \item java. Javascript version (at start of interview)
#'   \item age_cat. Year of birth category
#'   \item gender. Gender
#'   \item misconduct. Sum score of five binary items on student misconduct
#'   \item misconduct2. String of responses to five binary items on student misconduct
#'   \item field. Major field of study
#'   \item education. Type of study program
#'   \item semester. Current semester
#'   \item working. Working next to studying
#'   \item germanlang. German language skills
#'   \item riskattitude. Risk attitude (GSOEP 11-point scale)
#'   \item gpa. Current grade point average
#'   \item pressure. Studying is a lot of pressure
#'   \item stressed. Feeling very stressed in exams
#'   \item exams. Number of exams taken
#'   \item numberpapers. Number of papers handed in
#'   \item RRmodel. Randomized Response Model
#'   \item p1. Randomized Response parameter p1
#'   \item p2. Randomized Response parameter p2
#' }
#'
#' @docType data
#' @keywords datasets
#' @name ETHBE
#' @author Marc Hoeglinger, Ben Jann and Andreas Diekmann
#' @references \url{https://ideas.repec.org/p/bss/wpaper/8.html}
#' @usage data(ETHBE)
#' @format A data frame in long format with 21405 rows and 29 variables
NULL
