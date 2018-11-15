# Copyright (C) 2016-2018 by Xiaolei Lab
# 
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
# 
# http://www.apache.org/licenses/LICENSE-2.0
# 
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.



MVP.Version <- function(start=TRUE, width=60) {
##############################################################################################
# MVP: A Memory-efficient, Visualization-enhanced, and Parallel-accelerated Tool For GWAS
# Designed by Lilin Yin, Haohao Zhang, and Xiaolei Liu
# Reriten by Haohao Zhang
# Build date: Aug 30, 2017
# Last update: Oct 22, 2018
##############################################################################################
    
    welcome <- "Welcome to MVP"
    title   <- "A Memory-efficient, Visualization-enhanced, and Parallel-accelerated Tool For GWAS"
    authors <- "Authors: Lilin Yin, Haohao Zhang, and Xiaolei Liu"
    contact <- "Contact: xiaoleiliu@mail.hzau.edu.cn"
    logo_l  <- c(" __  ____      _______",
                 "|  \\/  \\ \\    / /  __ \\",
                 "| \\  / |\\ \\  / /| |__) |",
                 "| |\\/| | \\ \\/ / |  ___/",
                 "| |  | |  \\  /  | |",
                 "|_|  |_|   \\/   |_|")
    logo_s  <- c(" __  __  __   __  ___",
                 "|  \\/  | \\ \\ / / | _ \\",
                 "| |\\/| |  \\ V /  |  _/",
                 "|_|  |_|   \\_/   |_|")

    if (start) {
        print.info(welcome = welcome, title = title, logo = logo_s, authors = authors, contact = contact, line = '=', width = width)
    } else {
        rule(center = "MVP ACCOMPLISHED", width = width)
    }
}


#' Print R Package information, include title, short_title, logo, version, authors, contact
#'
#' @param welcome "Welcom to <Packagename>"
#' @param title long text to introduct package
#' @param short_title short label, top-left of logo
#' @param logo logo
#' @param version short label, bottom-right of logo
#' @param authors 
#' @param contact email or website
#' @param line 1, 2, or char
#' @param width banner width
#'
#' @author Haohao Zhang
#' @Build_date: Oct 22, 2018
#' @Last_update: Oct 22, 2018
print.info <- function(welcome=NULL, title=NULL, short_title=NULL, logo=NULL, version=NULL, authors=NULL, contact=NULL, line = 1, width=NULL) {
    msg <- c()
    # width
    if (is.null(width)) { width <- getOption('width') }
    # version
    if (is.null(version)) {
        if (getPackageName() == ".GlobalEnv") {
            version <- "devel"
        } else {
            version <- as.character(packageVersion(getPackageName()))
        }
    }
    # welcome
    if (is.null(welcome)) { 
        if (getPackageName() == ".GlobalEnv") {
            welcome <- ""
        } else {
            welcome <- paste0("Welcome to ", getPackageName())
        }
    }
    msg <- c(msg, rule(center = welcome, line = line, width = width))
    # title
    if (!is.null(title)) {
        msg <- c(msg, rule_wrap(string = title, width = width, align = "center"))
    }
    
    # align logo
    logo_width <- max(sapply(logo, nchar))
    for (i in 1:length(logo)) {
        # l <- stringr::str_pad(logo[i], "right", width = logo_width)
        l <- paste0(logo[i], paste(rep(" ", logo_width - nchar(logo[i])), collapse = ""))
        # blank <- stringr::str_match(logo[i], "^\\s+")
        l <- rule(center = l, line = " ", width = width)
        # if (!is.na(blank)) { l <- paste0(blank, substr(l, 1, nchar(l) - nchar(blank))) }
        msg <- c(msg, l)
    }
    
    # paste short_title label to logo top-left
    if (!is.null(short_title)) {
        i <- length(msg) - length(logo) + 1
        msg[i] <- paste_label(msg[i], paste0(short_title), side = "left")
    }
    
    # paste version label to logo bottom-right
    msg[length(msg)] <- paste_label(msg[length(msg)], paste0("Version: ", version), side = "right")
    
    # authors
    if (!is.null(authors)) {
        msg <- c(msg, rule_wrap(string = authors, align = "left", line = " ", width = width))
    }
    # contact
    if (!is.null(contact)) {
        msg <- c(msg, rule_wrap(string = contact, align = "left", line = " ", width = width))
    }
    # bottom line
    msg <- c(msg, rule(line = line, width = width))
    
    cat(msg, sep = "\n")
}


#' wrap text to multiple line, align left, right or center.
#' 
#' by using base::strwarp and cli::rule.
#' @author Haohao Zhang
#' @Build_date: Oct 22, 2018
#' @Last_update: Oct 22, 2018
rule_wrap <- function(string, width, align = "center", line = " ", margin = 2, ...) {
    # define
    msg <- c()
    lines <- strwrap(string, width = width - margin * 2)
    args <- list(width = width, line = line, ...)    # to cli::rule()
    
    # wrap
    for (i in 1:length(lines)) {
        args[[align]] <- lines[i]
        l <- do.call(rule, args)
        if (line == " " && align == "left" && substr(l, 1, 1) == " ") {
            l <- paste0(substring(l, 2), " ")
        }
        msg <- c(msg, l)
    }
    return(msg)
}


#' Paste label to a line
#' @author Haohao Zhang
#' @Build_date: Oct 22, 2018
#' @Last_update: Oct 22, 2018
paste_label <- function(line, label, side = "right", margin = 2) {
    if (side == "right") {
        end   <- nchar(line) - margin
        start <- end - (nchar(label) - 1)
    } else {
        start <- 1 + margin
        end   <- start + (nchar(label) - 1)
    }
    substr(line, start, end) <- label
    return(line)
}


format_time <- function(x)
{
    h <- x %/% 3600
    m <- (x %% 3600) %/% 60
    s <- ((x %% 3600) %% 60)
    index <- which(c(h, m, s) != 0)
    num <- c(h, m, s)[index]
    num <- round(num, 0)
    char <- c("h", "m", "s")[index]
    return(paste0(num, char, collapse = ""))
}