# Procedure to get bold_qc data:
library(httr2)

# Step 1: Query BOLD Systems API for Quebec records
query_request <- request("https://portal.boldsystems.org/api/query") |>
  req_url_query(
    query = "geo:province/state:Quebec",
    extent = "full"
  ) |>
  req_perform()

# Step 2: Download the data in TSV format
download_request <- request(glue::glue("https://portal.boldsystems.org/api/documents/{resp_body_json(query_request)$query_id}/download")) |>
  req_url_query(format = "tsv") |>
  req_perform()

# Extract the response body and write to file
bold_data <- download_request |> resp_body_string()
writeLines(bold_data, "results/bold_qc_data.tsv")


