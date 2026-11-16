MF.render_all(RES, META)
print(sorted(p.name for p in (RES / "figures").glob("Fig*.png")))