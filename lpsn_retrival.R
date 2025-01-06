
library(LPSN)



# Create an API access object

lpsn_access <- open_lpsn(username = "kamalavarshinibalakrishnan@gmail.com", password = "dsmz2001")

# Fetch details about Escherichia coli

prokaryote_info <- fetch(lpsn_access, ids = "Escherichia coli")
# Print the result
print(prokaryote_info)

# Convert to a data frame if needed
prokaryote_df <- as.data.frame(prokaryote_info)

# View the first few rows of the data frame
head(prokaryote_df)

