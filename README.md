gcloud builds submit --tag gcr.io/stavia-469704/stavia --project=stavia-469704

gcloud run deploy stavia --image gcr.io/stavia-469704/stavia --platform managed --project=stavia-469704 --allow-unauthenticated --memory 8Gi --region asia-east1 --timeout 3600 --cpu 2

gcloud run services list --project=stavia-469704

gcloud logs read --project=stavia-469704 --limit=50

gcloud run services logs read stavia --project=stavia-469704 --limit=50

gcloud run services describe stavia --project=stavia-469704

gcloud run services add-iam-policy-binding stavia --region asia-east1 --member=serviceAccount:414677286488-compute@developer.gserviceaccount.com --role=roles/storage.objectAdmin
