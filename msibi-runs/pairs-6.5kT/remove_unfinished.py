import signac

project = signac.get_project()
for job in project:
        if job.doc.done is False:
                job.remove()

