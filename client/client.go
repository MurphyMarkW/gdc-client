package client

import "io"
import "fmt"
import "net/http"

const RESOURCE = "https://%s:%d%s"

type Client struct {
	http.Client
	host string
	port int
}

func NewClient(host string, port int) *Client {
	return &Client{
		host: host,
		port: port,
	}
}

func (c *Client) Request(verb string, path string, token string, body io.Reader) (*http.Request, error) {
	var url = fmt.Sprintf("https://%s:%d%s", c.host, c.port, path)
	var req, err = http.NewRequest(verb, url, body)

	req.Header.Set("X-Auth-Token", token)

	return req, err
}

func (c *Client) Head(path string, token string) (*http.Response, error) {
	// TODO handle error cases
	var req, _ = c.Request(http.MethodHead, path, token, nil)

	return c.Do(req)
}

func (c *Client) Get(path string, token string) (*http.Response, error) {
	// TODO handle error cases
	var req, _ = c.Request(http.MethodGet, path, token, nil)

	return c.Do(req)
}

func (c *Client) Put(path string, token string, body io.Reader) (*http.Response, error) {
	// TODO handle error cases
	var req, _ = c.Request(http.MethodPut, path, token, body)

	return c.Do(req)
}

func (c *Client) Post(path string, token string, body io.Reader) (*http.Response, error) {
	// TODO handle error cases
	var req, _ = c.Request(http.MethodPost, path, token, body)

	return c.Do(req)
}

func (c *Client) Delete(path string, token string) (*http.Response, error) {
	// TODO handle error cases
	var req, _ = c.Request(http.MethodDelete, path, token, nil)

	return c.Do(req)
}

func (c *Client) Download(id string, token string) (io.ReadCloser, error) {
	const RESOURCE = "/v0/data/%s"
	var resource = fmt.Sprintf(RESOURCE, id)
	var resp, err = c.Get(resource, token)

	return resp.Body, err
}
