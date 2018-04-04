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

for batch in `cat batches`; do
    echo "at batch ${batch}"
    mkdir -p jobs/${batch}

    cat job_template.sh | sed "s/REPLACE_BATCH/${batch}/g" > jobs/${batch}/job_script.sh
    chmod 755 jobs/${batch}/job_script.sh

    echo "cd jobs/${batch}" >> submit.sh
    
    # Cluster specific submission.
    # Change as needed.
    echo "llsubmit job_script.sh" >> submit.sh
    echo "qsub job_script.sh" >> submit.sh
    echo "cd -" >> submit.sh
    echo ""     >> submit.sh
done
