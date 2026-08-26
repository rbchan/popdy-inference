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
  content <- gsub("<<(.*?)>>=", "```\\{r, \\1\\}", content)
  # Replace isolated '@' on their own line with standard markdown backticks
  content <- gsub("(?m)^\\s*@\\s*\\$", "```", content, perl = TRUE)
  content <- gsub("@", "```", content)
  
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
