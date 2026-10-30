"""The layers of the reconstruction and their default switches (no heavy imports, so that the
command line starts quickly)."""
LAYERS = ("surface", "sphere", "open", "closed", "streamer", "gcs", "radio", "c2", "cor1", "cor2", "metis")
ON = dict(surface=True, sphere=True, open=True, closed=True, streamer=True, gcs=True,
          radio=False, c2=False, cor1=False, cor2=False, metis=False)
