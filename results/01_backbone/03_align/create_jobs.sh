#!/bin/bash

# Safety. Do not want to delete important stuff.
if [ -d "jobs/" ]; then
    echo "Directory jobs/ already exists."
    read -p "Continue? [Y/n] " -n 1 -r response
    echo
    if [[ ! $response = "y" ]] && [[ ! $response = "" ]]; then
        echo "Aborted."
        exit
    fi
fi

# Fresh up the directories.
rm -rf "submit.sh"
rm -rf "jobs/"
mkdir -p "jobs"

echo -e "#!/bin/bash\n" > submit.sh
chmod 755 submit.sh

for batch in `ls batches`; do
    mkdir -p jobs/${batch}

    queries=`awk -vORS=, '{ print "${QUERYDIR}/"$1 }' batches/${batch} | sed 's/,$/\n/'`

    #echo $batch
    #echo $queries

    cat job_template.sh | sed "s/REPLACE_BATCH/${batch}/g" | sed "s?REPLACE_QUERIES?${queries}?g" > jobs/${batch}/job_script.sh
    chmod 755 jobs/${batch}/job_script.sh

    echo "cd jobs/${batch}" >> submit.sh
   
   	# Submit depending on cluster environment. 
   	# Change as needed
    echo "llsubmit job_script.sh" >> submit.sh
    echo "qsub job_script.sh" >> submit.sh

    echo "cd -" >> submit.sh
    echo ""     >> submit.sh
done

