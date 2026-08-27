rnw2qmd <- function(input_file, output_file = NULL) {
  if (is.null(output_file)) {
    output_file <- sub("\\.[Rr]nw\\$", ".qmd", input_file)
  }
  
  ## Read the full .Rnw content
  lines <- readLines(input_file, warn = FALSE)
  content <- paste(lines, collapse = "\n")
  
  ## 1. Strip the LaTeX preamble up to \begin{document}
  ## Inject a basic Quarto Revealjs YAML header
  yaml_header <- "---\n## title: \"Converted Presentation\"\nformat:\n  revealjs:\n    theme: default\n    pagetitle: 'Dog'\n    slide-number: true\nlang: en\n---\n\n"
  if (grepl("\\\\begin\\{document\\}", content)) {
    content <- sub("(?s).*?\\\\begin\\{document\\}", yaml_header, content, perl = TRUE)
  } else {
    content <- paste0(yaml_header, content)
  }
  
  ## 2. Remove standard closing LaTeX document tags
  content <- gsub("\\\\end\\{document\\}", "", content)
  
  ## 3. Convert \begin{frame}{Slide Title} to markdown header '## Slide Title'
  content <- gsub("\\\\begin\\{frame\\}\\[(.*?)\\]", "", content)
  ## content <- gsub("\\\\begin\\{frame\\}", "## Slide", content) # Fallback for title-less frames
  content <- gsub("\\\\begin\\{frame\\}", "", content) # Fallback for title-less frames
  content <- gsub("\\\\frametitle\\{(.*?)\\}", "## \\1", content) # Fallback for title-less frames
  content <- gsub("\\\\end\\{frame\\}", "", content)
  
  ## 4. Convert Sweave code chunks <<options>>= ... @ to ```{r, options} ... ```
  ## content <- gsub("<<(.*?)>>=", "```\\{r, \\1\\}", content)
  ## # Replace isolated '@' on their own line with standard markdown backticks
  ## content <- gsub("(?m)^\\s*@\\s*\\$", "```", content, perl = TRUE)
  ## content <- gsub("@", "```", content)

  # 4. Parse and translate Sweave chunks to Quarto chunks
  lines <- unlist(strsplit(content, "\n", fixed = TRUE))
  in_chunk <- FALSE
  new_lines <- c()
  
  for (line in lines) {
    # Match Sweave block start: <<label, option=value>>=
    ## if (!in_chunk && grepl("^\\s*<<", line) && grepl(">>=\\s*$", line)) {
    if (!in_chunk && grepl("<<", line) && grepl(">>", line)) {
        ## browser()
      in_chunk <- TRUE
      
      # Extract text inside << >>
      ## match_header <- regexec("^\\s*<<(.*)>>=\\s*$", line)
      match_header <- regexec("<<(.*)>>", line)
      header_content <- regmatches(line, match_header)[[1]]
      
      # Split by commas and trim whitespace
      opts <- unlist(strsplit(header_content[2], ",", fixed = TRUE))
      opts <- trimws(opts)
      
      label <- NULL
      quarto_opts <- c()
      
      for (i in seq_along(opts)) {
        opt <- opts[i]
        if (!grepl("=", opt, fixed = TRUE)) {
          if (i == 1) label <- opt # First nameless option is the label
        } else {
          kv <- unlist(strsplit(opt, "=", fixed = TRUE))
          key <- trimws(kv[1])
          val <- trimws(kv[2])
          # Convert booleans to lowercase for YAML compatibility
          if (tolower(val) == "true") val <- "true"
          if (tolower(val) == "false") val <- "false"
          quarto_opts <- c(quarto_opts, paste0("#| ", key, ": ", val))
        }
      }
      
      new_lines <- c(new_lines, "```{r}")
      if (!is.null(label) && label != "") {
        new_lines <- c(new_lines, paste0("#| label: ", label))
      }
      if (length(quarto_opts) > 0) {
        new_lines <- c(new_lines, quarto_opts)
      }
      
    ## } else if (in_chunk && grepl("^\\s*@\\s*\\$", line)) {
    } else if (in_chunk && grepl("@", line)) {
      # End of Sweave block: @
      in_chunk <- FALSE
      new_lines <- c(new_lines, "```")
    } else {
      # Document content or code body
      new_lines <- c(new_lines, line)
    }
  }
  
  content <- paste(new_lines, collapse = "\n")
  
  
  ## 5. Clean up standard list environments
  content <- gsub("\\\\begin\\{itemize\\}", "", content)
  content <- gsub("\\\\end\\{itemize\\}", "", content)
  content <- gsub("\\\\begin\\{enumerate\\}", "", content)
  content <- gsub("\\\\end\\{enumerate\\}", "", content)
  content <- gsub("\\\\item\\s*", "* ", content)

  ## 6. Clean up sections
  content <- gsub("\\\\section\\{(.*?)\\}", "# \\1", content) # Fallback for title-less frames

  ## 7. Columns
  content <- gsub("\\\\begin\\{columns\\}", "::: columns", content)
  content <- gsub("\\\\end\\{columns\\}", ":::", content)
  content <- gsub("\\\\begin\\{column\\}\\{(.*?)\\}", "::: \\{.column width='50\\%'\\}", content) 
  content <- gsub("\\\\end\\{column\\}", ":::", content)
  
  ## 8. Graphics
  content <- gsub("\\\\fbox\\{", "", content)
  content <- gsub("\\\\includegraphics", "", content)
#  content <- gsub("\\[(.*?)\\]", "\\[\\]", content)
#  content <- gsub("\\[\\]\\{(.*?)\\}", "\\!\\[\\](\\1)", content)

  ## 9. Misc
  content <- gsub("\\\\vfill", "", content)
  content <- gsub("\\\\hfill", "", content)
  content <- gsub("\\\\pause", "", content)
  content <- gsub("\\\\huge", "", content)
  content <- gsub("\\\\LARGE", "", content)
  content <- gsub("\\\\Large", "", content)
  content <- gsub("\\\\large", "", content)
  content <- gsub("\\\\normalsize", "", content)
  content <- gsub("\\\\footnotesizesize", "", content)
  content <- gsub("\\\\small", "", content)
  content <- gsub("\\\\centering", "", content)
  content <- gsub("\\\\begin\\{center\\}", "<center>", content)
  content <- gsub("\\\\end\\{center\\}", "</center>", content)
  content <- gsub("\\\\begin\\{align*\\}", "$$", content)
  content <- gsub("\\\\end\\{align*\\}", "$$", content)
  content <- gsub("\\\\vspace\\{(.*?)\\}", "", content) 
  content <- gsub("\\\\textbf\\{(.*?)\\}", "** \\1 **", content) 
  content <- gsub("\\\\textit\\{(.*?)\\}", "* \\1 *", content) 
  content <- gsub("\\\\tableofcontents", "", content)
  content <- gsub("^\\[\\]", "", content, perl=TRUE)
  content <- gsub("^\\s*\\[\\s*\\]\\s*$", "", content, perl=TRUE)
  content <- gsub("\\\\", "<br>", content, fixed=TRUE)
  
  # Write out the new .qmd file
  writeLines(content, output_file)
  message("Success! Converted file saved to: ", output_file)
}


### TODO: inr, ```r argument```
