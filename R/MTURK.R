#' MTurk Survey on "Mood and Personality"
#'
#' Data from an online validation experiment in which respondents' self-reports of norm breaking behavior were validated against observed actual behavior.
#' After playing a dice game, respondents were asked whether they played honestly, using one of several randomly assigned sensitive question techniques.
#' Furthermore, three other sensitive questions on shoplifting, tax evasion, and voting were asked.
#' The dataset contains the randomized responses from 6152 Amazon Mechanical Turk (MTURK) workers.
#' Each row holds the response to one question for one respondent.
#' The variables are as follows:
#'
#' \itemize{
#'   \item id. Identification code of the respondent
#'   \item Question. Which question was asked
#'   \item RR_response. Binary randomized response
#'   \item RRp1. Randomized Response parameter p1
#'   \item RRp2. Randomized Response parameter p2
#'   \item RRmodel. Randomized Response Model
#'   \item dicegame. Dice game assignment (1: prediction, 2: roll-a-six)
#'   \item cheater. The respondent is classified as honest or cheater if dice game assignment was 'roll-a-six'
#'   \item agecategory. Age category
#'   \item education. Level of education
#'   \item employment. Employment status
#'   \item locationinterview. Interview location
#'   \item extraversion. Extraversion score on a scale of 2-10
#'   \item agreeableness. Agreeableness score on a scale of 2-10
#'   \item conscientiousness. Conscientiousness score on a scale of 2-10
#'   \item neuroticism. Neuroticism score on a scale of 2-10
#'   \item openness. Openness score on a scale of 2-10
#'   \item gender. Gender (0: female, 1: male)
#'   \item age. Age in years
#'   \item privacyquestion1. How well are respondents' anonymity and privacy protected? (1: very poorly, 2: rather poorly, 3: moderately, 4: rather well, 5: very well)
#'   \item privacyquestion2. How likely could respondents' sensitive behavior be disclosed by this survey? (1: impossible, 2: not likely, 3: somewhat likely, 4: quite likely, 5: very likely)
#'   \item privacyquestion3. Does the special technique absolutely protect your answers? (1: not at all, 2: a little, 3: moderately, 4: quite a bit, 5: definitely)
#'   \item privacyquestion4. Do you think you properly followed the instructions for the special technique? (1: not at all, 2: a little, 3: moderately, 4: quite a bit, 5: definitely)
#'   \item privacyquestion5. Did you understand how the technique protects respondents? (1: not at all, 2: a little, 3: moderately, 4: quite a bit, 5: definitely)
#'   \item region. Region code
#'   \item country. Country
#' }
#'
#' @docType data
#' @keywords datasets
#' @name MTURK
#' @author Marc Hoeglinger and Ben Jann
#' @references \url{https://ideas.repec.org/p/bss/wpaper/17.html}
#' @usage data(MTURK)
#' @format A data frame in long format with 24594 rows and 26 variables
NULL
