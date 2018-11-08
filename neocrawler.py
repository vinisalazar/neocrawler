import scrapy


class EntrySpider(scrapy.Spider):
    name="entries"
    start_urls = [
        'http://www.bacterio.net/abiotrophia.html'
        # 'http://www.bacterio.net/-allnamesac.html',
        # 'http://www.bacterio.net/-allnamesdl.html',
        # 'http://www.bacterio.net/-allnamesmr.html',
        # 'http://www.bacterio.net/-allnamessz.html'
    ]

    def parse(self, response):
        for p in response.css('a[href*="https://www.ncbi.nlm.nih.gov/nuccore/"]::text'):
            yield {
                # 'genus': p.css('span.taxon-subhead-16s::text').extract_first()
                'genus': p.extract()
            }
