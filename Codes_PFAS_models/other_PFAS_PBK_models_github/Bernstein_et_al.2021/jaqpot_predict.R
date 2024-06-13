# Create a prediction instance
df = data.frame("admin.type" = "iv",
                "admin.dose" = 10, 
                "admin.time" = 0,
                "BW"=0.25, "BW.times" = 0,
                "F_unabs" = 0, "sex" = "M",  check.names = FALSE)
# Set the model ID of the uploaded Jaqpot model
modelID <- "B0TYi2AA7vC5DM7qfLV0"
# Log in the Jaqpot server
jaqpotr::login.cred()
# Acquire model predictions
predictions <- jaqpotr::jaqpot.predict( df = df, modelID = modelID)
# Plot an output
plot(predictions$predictions$time,predictions$predictions$Blood_total)