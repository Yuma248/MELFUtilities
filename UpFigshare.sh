#!/bin/bash
# This script is based on the example provided by figshare https://docs.figshare.com/ with practical modificiation. 

set -e
while [[ $# -gt 0 ]]; do
    key="$1"
    case $key in
        -i|--input)
        FILE_PATH="$2"
        shift
        shift
        ;;
        -u|--base_url)
        BASE_URL="$2"
        shift
        shift
        ;;
        -t|--token)
        ACCESS_TOKEN="$2"
        shift
        shift
        ;;
        -p|--project_id)
        PROJECT_ID="$2"
        shift
        shift
        ;;
        -f|--item_id)
        ITEM_ID="$2"
        shift
        shift
        ;;
       -n|--item_name)
        ITEM_NAME="$2"
        shift
        shift
        ;;
        *)
        shift
        ;;
    esac
done

if [ -z "$FILE_PATH" ] || [ -z "ACCESS_TOKEN" ]; then
        echo "This script will upload a file to a figshare repository, it needs the whole path of the file to be uploaded, a URL from figshare and a token to access the repository or project. If the file is not require to be uploaded in a specific project or item, you can use the general URL https://api.figshare.com/v2/account/articles. But if you need to put in specific project or item (folder) you will require the project ID (using the -p option) and/or the item ID (using the -f option). If you already upload one file in the project an item, you can integrate the item ID in the URL and add just /files at the end, see example"   
        echo -e "Usage:\nUpFigshare.sh\n\t-i <path to input file>\n\t-u <URL of figshare account/project/item, requiered if you have specific project and item ready, default https://api.figshare.com/v2/account/articles>\n\t-t <Access token for the figsahre account>\n\t-p <project ID, this is required just if you require to upload in specific project>\n\t-f <item ID, this is required just if you need to upload in specific item>\n\t-n <item name, this is required is you are creating a new item, default NEW UPLOAD>\n\n"
        echo -e "For example:\nUpFigshare.sh -i /sanpper/genome/genomev2.fasta -u https://api.figshare.com/v2/account/articles/22707367/files -t 75050303821b87ab7c72038ab9eaf02d853766bd8f7cc695f390d5b9cdeda1fd230c462a7cb7c7a67ef7c507f27f3fde647c41451574d54609ef477874c4z11\n\n"
        exit 1

fi

if [ -z "$BASE_URL" ]; then
	BASE_URL="https://api.figshare.com/v2/account/articles"
fi

if [ -z "$ITEM_NAME" ]; then
        ITEM_NAME="NEW UPLOAD"
fi


FILE_SIZE=$(stat -c%s $FILE_PATH)
MD5=($(md5sum $FILE_PATH))

FILE_NAME=$(echo $FILE_PATH | awk -F'/' '{print $NF}')
##echo "This file name is $FILE_NAME"


if echo "$BASE_URL" | grep -P -q -v  "files$"; then
	if [ -z "$ITEM_ID" ]; then 
		# Create a new item
		echo 'Creating a new item...'
		if [ -z "$PROJECT_ID" ] ; then
			RESPONSE=$(curl -s -f -d "{\"title\": \"$ITEM_NAME\"}" -H 'Authorization: token '$ACCESS_TOKEN -H 'Content-Type: application/json' -X POST "$BASE_URL") || { echo "Error: Creating new item failed, check that your URL and your token are correct"; exit 1; }
		else
			BASE_URL="https://api.figshare.com/v2/account/projects/${PROJECT_ID}/articles"
			RESPONSE=$(curl -s -f -d "{\"title\": \"$ITEM_NAME\"}" -H 'Authorization: token '$ACCESS_TOKEN -H 'Content-Type: application/json' -X POST "$BASE_URL") || { echo "Error: Creating new item failed, check that your URL, your token and your project ID are correct"; exit 1; }
			#echo "The location of the created item is "$RESPONSE
		fi
		echo ''
		# Retrieve item id
		echo 'Retrieving the item id...'
		ITEM_ID=$(echo "$RESPONSE" | sed -r "s/.*\/([0-9]+).*/\1/")
		echo "The item id is "$ITEM_ID" You can used this number to upload other files in the same folder"
		echo ''
	fi
	BASE_URL="https://api.figshare.com/v2/account/articles/${ITEM_ID}/files"
fi


# Initiate new upload:
echo 'A new upload had been initiated...'
RESPONSE=$(curl -s -f -d '{"md5": "'${MD5}'", "name": "'${FILE_NAME}'", "size": '${FILE_SIZE}'}' -H 'Content-Type: application/json' -H 'Authorization: token '$ACCESS_TOKEN -X POST "$BASE_URL") || { echo "Error: Creating new file failed, check that your URL and token are correct"; exit 1; }
#echo $RESPONSE
echo ''

# Retrieve file id
echo 'The file id is retrieved...'
FILE_ID=$(echo "$RESPONSE" | sed -r "s/.*\/([0-9]+).*/\1/")
#echo 'The file id is: '$FILE_ID
echo ''

# Retrieve the upload url
echo 'Retrieving the upload URL...'
RESPONSE=$(curl -s -f -H 'Authorization: token '$ACCESS_TOKEN -X GET "$BASE_URL/$FILE_ID")
UPLOAD_URL=$(echo "$RESPONSE" | sed -r 's/.*"upload_url":\s"([^"]+)".*/\1/')
#echo 'The upload URL is: '$UPLOAD_URL
echo ''

# Retrieve the upload parts
echo 'Retrieving the part value...'
RESPONSE=$(curl -s -f -H 'Authorization: token '$ACCESS_TOKEN -X GET "$UPLOAD_URL")
PARTS_SIZE=$(echo "$RESPONSE" | sed -r 's/"endOffset":([0-9]+).*/\1/' | sed -r 's/.*,([0-9]+)/\1/')
PARTS_SIZE=$(($PARTS_SIZE+1))
#echo 'The part value is: '$PARTS_SIZE
echo ''


# Split item into needed parts
echo 'Spliting the provided item into parts process had begun...'
split -b$PARTS_SIZE $FILE_PATH part_ --numeric=1

#echo 'Process completed!'

# Retrive the number of parts
MAX_PART=$((($FILE_SIZE+$PARTS_SIZE-1)/$PARTS_SIZE))
#echo 'The number of parts is: '$MAX_PART
echo ''

# Perform the PUT operation of parts
echo 'Perform the PUT operation of parts...'
for ((i=1; i<=$MAX_PART; i++))
do
    PART_VALUE='part_'$i
    if [ "$i" -le 9 ]
    then
        PART_VALUE='part_0'$i
    fi
    RESPONSE=$(curl -s -f -H 'Authorization: token '$ACCESS_TOKEN -X PUT "$UPLOAD_URL/$i" --data-binary @$PART_VALUE)
    echo "Done uploading part nr: $i/"$MAX_PART
done

echo 'Process was finished!'
echo ''

# Complete upload
echo 'Completing the file upload...'
RESPONSE=$(curl -s -f -H 'Authorization: token '$ACCESS_TOKEN -X POST "$BASE_URL/$FILE_ID") || { echo "Error: Uploading file failed, check that your URL and token are correct"; exit 1; }
echo 'Done!'
echo ''

#remove the part files
rm part_*

# List all of the existing items
#RESPONSE=$(curl -s -f -H 'Authorization: token '$ACCESS_TOKEN -X GET "$BASE_URL")
#echo 'New list of items: '$RESPONSE
#echo ''

