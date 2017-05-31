from skbuild import setup
import subprocess

setup(
    name="invdet",
    version=("0.1."+subprocess.check_output(["git", "rev-list", "--count", "9e02a86bcbb0698f327ed694cf602e5e477174d3..HEAD"])).rstrip(),
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
