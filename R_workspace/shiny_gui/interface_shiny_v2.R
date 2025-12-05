# Install required packages if not already installed
if (!requireNamespace("shiny", quietly = TRUE)) install.packages("shiny")
if (!requireNamespace("rhandsontable", quietly = TRUE)) install.packages("rhandsontable")
if (!requireNamespace("DT", quietly = TRUE)) install.packages("DT")
if (!requireNamespace("ggplot2", quietly = TRUE)) install.packages("ggplot2")
if (!requireNamespace("ggpubr", quietly = TRUE)) install.packages("ggpubr")


library(shiny)
library(rhandsontable)
library(DT)
library(ggplot2)
library(ggpubr)

ui<-fluidPage(
  titlePanel("Prepare to Anger!"),
  sidebarLayout(
    sidebarPanel(
      actionButton("add_read", "Add Read Element"),
      actionButton("remove_read", "Remove Read Element"),
      actionButton("add_path", "Add Path"),
      actionButton("remove_path", "Remove Path"),
      actionButton("submit", "Submit")
    ),
    mainPanel(
      rHandsontableOutput("read_table", width = "100%", height = "100%"),
      DTOutput("path_table"),
      plotOutput("read_layout_plot")
    )
  )
)

server<-function(input, output, session){
  read_data <- reactiveVal(data.frame(id = integer(0), seq = character(0), expected_length = integer(0),
    type = factor(levels = c("barcode", "adapter", "extraneous", "UMI", "read", "TSO", "poly_a", "poly_t")),
    stringsAsFactors = FALSE))
  path_data <- reactiveVal(data.frame(path_type = character(0), path_value = character(0),
    stringsAsFactors = FALSE))
  observeEvent(input$path_table_cell_edit,{
    info <- input$path_table_cell_edit
    modified_data <- path_data()
    modified_data[info$row, info$col]<-info$value
    path_data(modified_data)
  })
  observeEvent(input$add_read,{
    new_row<-data.frame(id = NA, seq = NA, expected_length = NA, type = NA,
      stringsAsFactors = FALSE)
    read_data(rbind(read_data(), new_row))
  })
  observeEvent(input$remove_read,{
    data <- read_data()
    if (nrow(data) > 0) {
      read_data(data[-nrow(data), , drop = FALSE])
    }
  })
  observeEvent(input$add_path,{
    new_row <- data.frame(path_type = NA, path_value = NA, stringsAsFactors = FALSE)
    current_data <- path_data()
    updated_data <- rbind(current_data, new_row)
    path_data(updated_data)
  })
  observeEvent(input$remove_path,{
    current_data<-path_data()
    if (nrow(current_data) > 0) {
      updated_data<-current_data[-nrow(current_data), , drop = FALSE]
      path_data(updated_data)
    }
  })
  output$path_table<-renderDT({
    datatable(path_data(), editable = TRUE)
  })
  observeEvent(input$submit,{
    # Display notification
    showNotification(paste0("Read Layout Data: ", toString(read_data())))
    showNotification(paste0("Paths Data: ", toString(path_data())))

    # Assign data to global variables
    assign("read_layout", read_data(), envir = .GlobalEnv)
    assign("path_layout", path_data(), envir = .GlobalEnv)
  })

  output$read_table<-renderRHandsontable({
    rhandsontable(read_data(), rowHeaders = NULL, width = 800, height = 300) %>%
      hot_col("id", type = "text") %>%
      hot_col("seq", type = "text") %>%
      hot_col("expected_length", type = "numeric") %>%
      hot_col("type", type = "dropdown", source = c("barcode", "adapter", "extraneous", "UMI", "read", "TSO", "poly_a", "poly_t"))
  })

  observe({
    if (!is.null(input$read_table)){
      read_data(hot_to_r(input$read_table))
    }
  })

  my_ggplot<-reactive({
    data <- read_data()
    if (nrow(data) > 0) {
      if (any(data$type == "read", na.rm = TRUE)) {
        data$expected_length[data$type == "read"] <- max(data$expected_length, na.rm = TRUE)
      }
      # Create variables for the x-axis (position in the read layout)
      data$xmin <- cumsum(c(0, head(data$expected_length, -1)))
      data$xmax <- cumsum(data$expected_length)
      # Manually specify colors for each type
      color_map<-c("barcode" = "#70AD47", "adapter" = "black", "extraneous" = "gold", "UMI" = "orange", "read" = "lightgrey", "TSO" = "#20B7AD","poly_a" = "red", "poly_t" = "#008BD3")
      plot_obj<-ggplot2::ggplot(data, aes(xmin = xmin, xmax = xmax, ymin = 0.4, ymax = 0.6, fill = type)) +
        geom_rect(na.rm = TRUE)+
        scale_fill_manual(values = color_map)+
        theme_void()
      plot_obj<-ggpar(p = plot_obj,
        main = "Read Layout", submain = "5' to 3'",
        font.main = c("36", "bold"), font.submain = c(24,"bold"),
        ylim = c(0,1), xlab = FALSE, ylab = FALSE, ticks = FALSE, tickslab = FALSE)
      plot_obj<-plot_obj+theme(plot.title = element_text(hjust = 0.5), plot.subtitle = element_text(hjust = 0.5))
      return(plot_obj)
    }
  })
  output$read_layout_plot <- renderPlot({
    my_ggplot()
  })
}

shinyApp(ui, server)
