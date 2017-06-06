from skbuild import setup
import subprocess

setup(
    name="invdet",
    version=("0.3."+subprocess.check_output(["git", "rev-list", "--count", "59926c99ea38e4eccd7e44912d99a275556d7a2b..HEAD"])).rstrip(),
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
