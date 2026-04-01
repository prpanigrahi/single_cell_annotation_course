# Install dependencies
#install.packages("devtools")
# install.packages("reticulate")

# Install CASSIA
# devtools::install_github("ElliotXie/CASSIA/CASSIA_R")

library(reticulate)
library(CASSIA)

#
my_env_path <- "./venv"
virtualenv_create(envname = my_env_path)
use_virtualenv(my_env_path, required = TRUE)
setup_cassia_env()

# For local LLMs - no API key needed (e.g., Ollama)
setLLMApiKey(provider = "http://localhost:11434/v1", persist = TRUE)

top10 = top_markers %>% group_by(cluster) %>% slice_max(order_by = avg_log2FC, n = 10)

# Core annotation
runCASSIA_batch(
  marker = top10,                # Marker data from FindAllMarkers
  output_name = "cassia_results",              # Output file name
  tissue = "Blood",                  # Tissue type
  species = "Human",                           # Species
  provider = "http://localhost:11434/v1",
  model = "gpt-oss:20b",                     # API provider
  max_workers = 4                              # Number of parallel workers
)



