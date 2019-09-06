package download

import client "github.com/MurphyMarkW/gdc-client/client"

import "fmt"
import "net/http"

const AGENT = "go-gdc/0.0.0 (%s) %s (%s)"
const RESOURCE = "https://%s:%d%s"

type DownloadClient struct {
	client.Client
}

func NewDownloadClient(host string, port uint16) DownloadClient {
	return DownloadClient{*client.NewClient(host, port)}
}

func (c *DownloadClient) Download(id string, token string) (*http.Response, error) {
	const RESOURCE = "/v0/data/%s"
	var resource = fmt.Sprintf(RESOURCE, id)
	var res, err = c.Get(resource, token)

	return res, err
}
