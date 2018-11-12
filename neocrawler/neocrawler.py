import scrapy
from urllib.parse import urljoin

"""
This file crawls the Bacterio 'names' pages. It should open the link for
each Genus page, and then extract the accession number for each one.
"""

class BacteriumSpider(scrapy.Spider):
    name="bacterium"

    start_urls = [
        'http://www.bacterio.net/-allnamesac.html',
        'http://www.bacterio.net/-allnamesdl.html',
        'http://www.bacterio.net/-allnamesmr.html',
        'http://www.bacterio.net/-allnamessz.html'
    ]

    def parse(self, response):
        bacteria = response.css('a[href$=".html"]').extract()
        bacteria = [b for b in bacteria if '¤' in b]
        for b in bacteria:
            b = b.split('"')[1]
            url = urljoin(response.url, b)
            yield scrapy.Request(url, callback=self.parse_bacteria)

    def parse_bacteria(self, response):
        for p in response.css(
            'a[href*="https://www.ncbi.nlm.nih.gov/"]::text'
            ):
            yield { 'ncbi_id': p.extract() }
