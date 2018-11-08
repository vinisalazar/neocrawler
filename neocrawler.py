import scrapy


class EntrySpider(scrapy.Spider):
    name="entries"
    start_urls = [
        'http://www.bacterio.net/-allnamesac.html'
    ]
    def parse(self, response):
        for p in response.css('div.main-text'):
            yield {
                'genus': p.css('span.class::genusspecies')
            }
