'''
Some extensions to brp for plotting single pulse diagnostics.
'''
import StringIO
from base64 import encodestring

from PIL import Image, ImageDraw

from brp.svg.et_import import ET
from brp.svg.plotters.scatter import ScatterPlotter
from brp.svg.colornames import svg_color2rgba_color

# -----------------------------------------------------------------------------
# -- Extensions to brp for plotting single pulse candidates -------------------


class DetectionSymbol(object):
    def __init__(self, *args, **kwargs):
        self.color = kwargs.get('color', 'black')
        self.link = kwargs.get('link', '')

    def draw_xy(self, root_element, x, y, **kwargs):
        color = kwargs.get('color', self.color)
        sigma = kwargs['sigma']
        radius = min(0.5 * (sigma - 5) + 1, 8)
        link = kwargs.get('link', self.link)

        if link:
            root_element = ET.SubElement(root_element, 'a')
            root_element.set('xlink:href', link)

        p = ET.SubElement(root_element, 'circle')
        p.set('cx', '%.2f' % x)
        p.set('cy', '%.2f' % y)
        p.set('r', '%.2f' % radius)
        p.set('stroke', color)
        p.set('fill', 'none')

    def draw(self, root_element, x_transform, y_transform, x, y, **kwargs):
        nx = x_transform(x)
        ny = y_transform(y)
        self.draw_xy(root_element, nx, ny, **kwargs)

    def rdraw(self, imdraw, x_transform, y_transform, x, y, **kwargs):
        nx = x_transform(x)
        ny = y_transform(y)

        rgba_color = svg_color2rgba_color(self.color)
        sigma = kwargs['sigma']
        radius = min(0.5 * (sigma - 5) + 1, 8)

        bbox = [
            nx - radius,
            ny - radius,
            nx + radius,
            ny + radius,
        ]
        imdraw.ellipse(bbox, fill=None, outline=rgba_color)


class DetectionPlotterTuples(ScatterPlotter):
    '''
    Expexts (dm, sigma, t, sample, downfact) tuples.

    Note also works with the new style:
    (dm, sigma, t, sample, downfact, t - delta_t, t + delta_t, i)
    '''
    def __init__(self, candidates, *args, **kwargs):
        self.candidates = candidates[:]
        self.color = kwargs.get('color', 'black')
        self.symbols = kwargs.get('symbols',
                                  [DetectionSymbol(color=self.color)])

    def prepare_bbox(self, data_bbox=None):
        if data_bbox is None:
            tmp = self.candidates[0]
            data_bbox = [tmp[2], tmp[0], tmp[2], tmp[0]]
        else:
            data_bbox = list(data_bbox)

        for c in self.candidates:
            if c[2] < data_bbox[0]:
                data_bbox[0] = c[2]
            elif c[2] > data_bbox[2]:
                data_bbox[2] = c[2]
            if c[0] < data_bbox[1]:
                data_bbox[1] = c[0]
            elif c[0] > data_bbox[3]:
                data_bbox[3] = c[0]
        return data_bbox

    def draw(self, root_element, x_transform, y_transform):
        for candidate in self.candidates:
            for symbol in self.symbols:
                symbol.draw(root_element, x_transform, y_transform,
                            candidate[2], candidate[0], sigma=candidate[1],
                            downfact=candidate[4])

    def rdraw(self, root_element, x_transform, y_transform, svg_bbox):

        width = svg_bbox[2] - svg_bbox[0]
        height = svg_bbox[1] - svg_bbox[3]
        assert width > 0
        assert height > 0

        im = Image.new('RGBA', (width, height), (255, 255, 255, 0))
        imdraw = ImageDraw.Draw(im)
        # ADD PLOT STUFF HERE
        for candidate in self.candidates:
            for symbol in self.symbols:
                symbol.rdraw(imdraw, x_transform, y_transform,
                             candidate[2], candidate[0], sigma=candidate[1],
                             downfact=candidate[4])
        # BOILERPLATE:
        tmp = StringIO.StringIO()
        im.save(tmp, format='png')
        encoded_png = 'data:image/png;base64,\n' + encodestring(tmp.getvalue())

        img = ET.SubElement(root_element, 'image')
        img.set('xlink:href', encoded_png)
        img.set('x', '%.2f' % min(svg_bbox[0], svg_bbox[2]))
        img.set('y', '%.2f' % min(svg_bbox[1], svg_bbox[3]))
        img.set('width', '%.2f' % width)
        img.set('height', '%.2f' % height)
        img.set('preserveAspectRatio', 'none')
