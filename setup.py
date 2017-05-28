from skbuild import setup
import subprocess

setup(
    name="invdet",
    version="0.2."+subprocess.check_output(["git", "rev-parse", "HEAD"])[:7],
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
