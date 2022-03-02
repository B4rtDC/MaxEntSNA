This repository contains the code associated with the paper "Maximum entropy networks for large scale social network node analysis". The scripts are written in Julia (>= 1.6) and Python (>=3.8).

To use the different scripts, you need to have downloaded the archives from the [Twitter information operations repository](https://transparency.twitter.com/en/reports/information-operations.html) 


The steps to undertake for each dataset in the Twitter information operations repository are the following:
1. extract the external tweets that need to be downloaded
2. download the external tweets with a tool of choice (e.g. [twarc](https://twarc-project.readthedocs.io/en/latest/), which we used here or [Hydrator](https://github.com/DocNow/hydrator))
3. build the relevant interaction network(s)
4. identify the significant interactions and project the interaction network on the user layer (using the [NEMtropy](https://github.com/nicoloval/NEMtropy) package)
5. analyze the results



## Obtaining the data
1. Download the files associated with a specific dataset. This will give you a set of files that holds the tweets and the user info that has the following naming conventions:
    ```julia
    dataset_users_csv_hashed.csv # holds the flagged user info
    dataset_tweets_csv_hashed_XX.csv # holds the tweets, can be multiple files indicated by the number XX
    ```
2. From the Analyse module, run the function `get_externals`. This will provide two files, one for the external retweets and one for the external replies to be downloaded.
3. Set your credentials for the Twitter API in `credentials.py`
4. Run the the function `harvestfiles` from `hydratorpipeline.py`. This will identify the files that hold the external data. The tweets will be downloaded into a .jsonl file. At the same time, a log file will be generated to track what percentage of messages could be recovered.

## Obtaining the networks
A Pluto notebook is provided that allows you to replicate the entire process. The results are also available in a webpage, should you not have the data.

## Exporting the network to Gephi
You can export each network for importation in Gephi by exporting the nodelist and edgelists.
