# R-shiny-app---Differential-expression-analysis

This app was developed to analyze the data of "mRNA-Seq Expression profiling of human post-mortem BA9 brain tissue for Huntington's Disease and neurologically normal individuals" which is publically available on the GEO browser at GSE64810

1. Libraries Used
The application utilizes several R libraries to handle data manipulation, visualization, and user interface enhancement.

Key libraries include:
shiny for creating the web app.
colourpicker and bslib for enhancing user interface elements.
DT, ggplot2, plotly, and similar packages for data visualization.
tidyverse for data manipulation.
Other specialized libraries like igraph, gplots, heatmaply, beeswarm, and vioplot for specific types of plots.

2. User Interface (UI)
   
The UI is structured into four main tabs, each designed for specific aspects of the data analysis:
Sample Information: Allows users to upload sample information through a CSV file and view summaries, detailed tables, and histograms of the sample data.
Counts Matrix: Focuses on uploading and exploring a counts matrix, including variance filtering, visualization of diagnostic plots (e.g., scatter plots of median variance and zeros), heatmap generation, and PCA plots.
Differential Expression: Supports uploading differential expression results for further analysis. Features include dynamic selection of axes for volcano plots, customizable color schemes, and interactive sliders for adjusting the p-value coloring threshold.
Individual Gene Expressions: Offers functionality for uploading normalized counts and sample information to generate plots (bar, box, beeswarm, violin) comparing gene expressions across different categories.

3. Server Logic
   
The server component implements the application's functionality, including:
Reading and processing uploaded files.
Generating dynamic summaries and plots based on user inputs and selected parameters.
Filtering data according to user-defined criteria (e.g., variance and nonzero samples threshold).
Creating interactive visualizations like PCA plots, volcano plots, and gene expression plots.

4. Features and Functionalities
   
File Uploads: Users can upload CSV files of sample information, counts matrix, and differential expression results.
Data Exploration: Tools for detailed exploration of the data include dynamic tables, summary statistics, histograms, and diagnostic plots.
Interactive Visualizations: Includes heatmap for gene expression, PCA plots for dimensionality reduction insights, volcano plots for differential expression analysis, and specific gene expression plots.
Customization: Users can customize plot appearances, select genes for analysis, and adjust thresholds for data filtering.
