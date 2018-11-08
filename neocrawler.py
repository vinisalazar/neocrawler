import scrapy


class EntrySpider(scrapy.Spider):
    name="entries"
    urls = [
        'http://www.bacterio.net/-allnamesac.html',
        'http://www.bacterio.net/-allnamesdl.html',
        'http://www.bacterio.net/-allnamesmr.html',
        'http://www.bacterio.net/-allnamessz.html'
    ]

    def parse(self, response):
        for p in response.css('p'):
            yield {
                'genus': p.css('span.class::genusspecies')
            }
