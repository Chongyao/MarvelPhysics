;;; Directory Local Variables
;;; For more information see (info "(emacs) Directory Variables")


((prog-mode . ((eval . (let ((root (projectile-project-root)))
                         (setq-local local-project-include-path
                                     (list
                                      (concat root "src/")
                                      (concat root "external/")
                                      (concat root "src/utils/")))
                         (setq-local company-c-headers-path-user
                                     (append company-c-headers-path-user
                                             local-project-include-path))
                         (setq-local flycheck-clang-include-path
                                     (append flycheck-clang-include-path
                                             local-project-include-path))
                         (defun concat_I (dirs)
                           "concat -I to the elements in the dirs"
                           (setq-local concat_res (list))
                           (while dirs
                             (setq-local concat_res (append concat_res (list (concat "-I" (car dirs)))))
                             (setq-local dirs (cdr dirs)))
                           concat_res)
                         (setq-local company-clang-arguments
                                     (append company-clang-arguments
                                             (concat_I local-project-include-path))))))))
