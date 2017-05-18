from skbuild import setup

setup(
    name="invdet",
    version="0.1.31",
    description="Invertion Detection package",
    author="zijuexiansheng",
    license="MIT",
    packages=["invdet", "invdet.util"],
    package_dir={"invdet": "source/pyinvdet",
                 "invdet.util": "source/pyinvdet/util"},
    entry_points={
            'console_scripts': ['invdet=invdet.invdet_main:main']
        }
)
