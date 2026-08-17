rnw2qmd <- function(input_file, output_file = NULL) {
  if (is.null(output_file)) {
    output_file <- sub("\\.[Rr]nw\\$", ".qmd", input_file)
  }
  
  ## Read the full .Rnw content
  lines <- readLines(input_file, warn = FALSE)
  content <- paste(lines, collapse = "\n")
  
  ## 1. Strip the LaTeX preamble up to \begin{document}
  ## Inject a basic Quarto Revealjs YAML header
  yaml_header <- "---\ntitle: \"Converted Presentation\"\nformat:\n  revealjs:\n    theme: default\n    slide-number: true\nlang: en\n---\n\n"
  if (grepl("\\\\begin\\{document\\}", content)) {
    content <- sub("(?s).*?\\\\begin\\{document\\}", yaml_header, content, perl = TRUE)
  } else {
    content <- paste0(yaml_header, content)
  }
  
  ## 2. Remove standard closing LaTeX document tags
  content <- gsub("\\\\end\\{document\\}", "", content)
  
  ## 3. Convert \begin{frame}{Slide Title} to markdown header '## Slide Title'
  ## content <- gsub("\\\\begin\\{frame\\}\\{(.*?)\\}", "## \\1", content)
  ## content <- gsub("\\\\begin\\{frame\\}", "## Slide", content) # Fallback for title-less frames
  content <- gsub("\\\\begin\\{frame\\}", "", content) # Fallback for title-less frames
  content <- gsub("\\\\frametitle\\{(.*?)\\}", "## \\1", content) # Fallback for title-less frames
  content <- gsub("\\\\end\\{frame\\}", "", content)
  
  ## 4. Convert Sweave code chunks <<options>>= ... @ to ```{r, options} ... ```
  content <- gsub("<<(.*?)>>=", "```\\{r, \\1\\}", content)
  # Replace isolated '@' on their own line with standard markdown backticks
  content <- gsub("(?m)^\\s*@\\s*\\$", "```", content, perl = TRUE)
  
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
  content <- gsub("\\[(.*?)\\]", "\\[\\]", content)
  content <- gsub("\\[\\]\\{(.*?)\\}", "\\!\\[\\](\\1)", content)
  
  
  # Write out the new .qmd file
  writeLines(content, output_file)
  message("Success! Converted file saved to: ", output_file)
}
