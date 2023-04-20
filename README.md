# enhancerGenePro
It takes a little longer to clone this repository as supporting files are of high lengths.
## Build the image 
sudo  docker build -t my-app-image .

## To run the server
sudo docker run -dp 8080:8080 my-app-image

## Server should be available on : localhost:8080

## To check logs
docker logs -f < containerid >

## To get container ID
docker ps
