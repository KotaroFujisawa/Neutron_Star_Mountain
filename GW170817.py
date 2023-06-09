import gwpy
from gwpy.timeseries import TimeSeries
template = TimeSeries.fetch_open_data("H1", "August 16, 2017", "August 17, 2017")
plot = template.plot()
ax = plot.gca()
ax.set_ylabel('Strain amplitude')
ax.set_title('GW170817 template')
plot.show()
